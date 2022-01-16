/**
 * K.jpg's OpenSimplex 2, faster variant
 *
 * - 2D is standard simplex implemented using a lookup table.
 * - 3D is "Re-oriented 4-point BCC noise" which constructs a
 *   congruent BCC lattice in a much different way than usual.
 * - 4D constructs the lattice as a union of five copies of its
 *   reciprocal. It successively finds the closest point on each.
 *
 * Multiple versions of each function are provided. See the
 * documentation above each, for more info.
 */
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "OpenSimplex2F.h"

#define PSIZE (2048)
#define PMASK (2047)
static const double N2 = 0.01001634121365712;
static const double N3 = 0.030485933181293584;
static const double N4 = 0.009202377986303158;

#define B0000 0x00
#define B0001 0x01
#define B0010 0x02
#define B0011 0x03
#define B0100 0x04
#define B0101 0x05
#define B0110 0x06
#define B0111 0x07
#define B1000 0x08
#define B1001 0x09
#define B1010 0x0a
#define B1011 0x0b
#define B1100 0x0c
#define B1101 0x0d
#define B1110 0x0e
#define B1111 0x0f

struct Grad2 {
	double dx, dy;
};

struct Grad3 {
	double dx, dy, dz;
};

struct Grad4 {
	double dx, dy, dz, dw;
};

struct LatticePoint2D {
	int xsv, ysv;
	double dx, dy;
};

static struct LatticePoint2D *new_LatticePoint2D(int xsv, int ysv)
{
	struct LatticePoint2D *this = calloc(1, sizeof(*this));

	this->xsv = xsv; this->ysv = ysv;
	double ssv = (xsv + ysv) * -0.211324865405187;
	this->dx = -xsv - ssv;
	this->dy = -ysv - ssv;
	return this;
}

struct LatticePoint3D {
	double dxr, dyr, dzr;
	int xrv, yrv, zrv;
	struct LatticePoint3D *nextOnFailure, *nextOnSuccess;
};

static struct LatticePoint3D *new_LatticePoint3D(int xrv, int yrv, int zrv, int lattice)
{
	struct LatticePoint3D *this = calloc(1, sizeof(*this));

	this->dxr = -xrv + lattice * 0.5; this->dyr = -yrv + lattice * 0.5; this->dzr = -zrv + lattice * 0.5;
	this->xrv = xrv + lattice * 1024; this->yrv = yrv + lattice * 1024; this->zrv = zrv + lattice * 1024;
	return this;
}

static void find_unique_pointers(struct LatticePoint3D *tree, struct LatticePoint3D *list[], int *n)
{
	if (!tree)
		return;
	int i, found;

	found = 0;
	for (i = 0; i < *n; i++) {
		if (list[i] == tree) {
			found = 1;
			break;
		}
	}
	if (!found) {
		list[*n] = tree;
		(*n)++;
	}
	find_unique_pointers(tree->nextOnFailure, list, n);
	find_unique_pointers(tree->nextOnSuccess, list, n);
}

static void free_LatticePoint3D(struct LatticePoint3D *list[], int n)
{
	int i;

	for (i = 0; i < n; i++)
		free(list[i]);
}

struct LatticePoint4D {
	int xsv, ysv, zsv, wsv;
	double dx, dy, dz, dw;
	double xsi, ysi, zsi, wsi;
	double ssiDelta;
};

static struct LatticePoint4D *new_LatticePoint4D(int xsv, int ysv, int zsv, int wsv)
{
	struct LatticePoint4D *this = calloc(1, sizeof(*this));

	this->xsv = xsv + 409; this->ysv = ysv + 409; this->zsv = zsv + 409; this->wsv = wsv + 409;
	double ssv = (xsv + ysv + zsv + wsv) * 0.309016994374947;
	this->dx = -xsv - ssv;
	this->dy = -ysv - ssv;
	this->dz = -zsv - ssv;
	this->dw = -wsv - ssv;
	this->xsi = 0.2 - xsv;
	this->ysi = 0.2 - ysv;
	this->zsi = 0.2 - zsv;
	this->wsi = 0.2 - wsv;
	this->ssiDelta = (0.8 - xsv - ysv - zsv - wsv) * 0.309016994374947;
	return this;
}

struct OpenSimplex2F_context {
	int16_t *perm;
	struct Grad2 *permGrad2;
	struct Grad3 *permGrad3;
	struct Grad4 *permGrad4;
};

#define ARRAYSIZE(x) (sizeof(x) / sizeof((x)[0]))

static struct Grad2 GRADIENTS_2D[PSIZE];
static struct Grad3 GRADIENTS_3D[PSIZE];
static struct Grad4 GRADIENTS_4D[PSIZE];
static struct LatticePoint2D *LOOKUP_2D[4];
static struct LatticePoint3D *LOOKUP_3D[8];
static struct LatticePoint4D *VERTICES_4D[16];

static struct Grad2 grad2[] = {
	{ 0.130526192220052,  0.99144486137381},
	{ 0.38268343236509,   0.923879532511287},
	{ 0.608761429008721,  0.793353340291235},
	{ 0.793353340291235,  0.608761429008721},
	{ 0.923879532511287,  0.38268343236509},
	{ 0.99144486137381,   0.130526192220051},
	{ 0.99144486137381,  -0.130526192220051},
	{ 0.923879532511287, -0.38268343236509},
	{ 0.793353340291235, -0.60876142900872},
	{ 0.608761429008721, -0.793353340291235},
	{ 0.38268343236509,  -0.923879532511287},
	{ 0.130526192220052, -0.99144486137381},
	{-0.130526192220052, -0.99144486137381},
	{-0.38268343236509,  -0.923879532511287},
	{-0.608761429008721, -0.793353340291235},
	{-0.793353340291235, -0.608761429008721},
	{-0.923879532511287, -0.38268343236509},
	{-0.99144486137381,  -0.130526192220052},
	{-0.99144486137381,   0.130526192220051},
	{-0.923879532511287,  0.38268343236509},
	{-0.793353340291235,  0.608761429008721},
	{-0.608761429008721,  0.793353340291235},
	{-0.38268343236509,   0.923879532511287},
	{-0.130526192220052,  0.99144486137381}
};

static struct Grad3 grad3[] = {
	{-2.22474487139,      -2.22474487139,      -1.0},
	{-2.22474487139,      -2.22474487139,       1.0},
	{-3.0862664687972017, -1.1721513422464978,  0.0},
	{-1.1721513422464978, -3.0862664687972017,  0.0},
	{-2.22474487139,      -1.0,                -2.22474487139},
	{-2.22474487139,       1.0,                -2.22474487139},
	{-1.1721513422464978,  0.0,                -3.0862664687972017},
	{-3.0862664687972017,  0.0,                -1.1721513422464978},
	{-2.22474487139,      -1.0,                 2.22474487139},
	{-2.22474487139,       1.0,                 2.22474487139},
	{-3.0862664687972017,  0.0,                 1.1721513422464978},
	{-1.1721513422464978,  0.0,                 3.0862664687972017},
	{-2.22474487139,       2.22474487139,      -1.0},
	{-2.22474487139,       2.22474487139,       1.0},
	{-1.1721513422464978,  3.0862664687972017,  0.0},
	{-3.0862664687972017,  1.1721513422464978,  0.0},
	{-1.0,                -2.22474487139,      -2.22474487139},
	{ 1.0,                -2.22474487139,      -2.22474487139},
	{ 0.0,                -3.0862664687972017, -1.1721513422464978},
	{ 0.0,                -1.1721513422464978, -3.0862664687972017},
	{-1.0,                -2.22474487139,       2.22474487139},
	{ 1.0,                -2.22474487139,       2.22474487139},
	{ 0.0,                -1.1721513422464978,  3.0862664687972017},
	{ 0.0,                -3.0862664687972017,  1.1721513422464978},
	{-1.0,                 2.22474487139,      -2.22474487139},
	{ 1.0,                 2.22474487139,      -2.22474487139},
	{ 0.0,                 1.1721513422464978, -3.0862664687972017},
	{ 0.0,                 3.0862664687972017, -1.1721513422464978},
	{-1.0,                 2.22474487139,       2.22474487139},
	{ 1.0,                 2.22474487139,       2.22474487139},
	{ 0.0,                 3.0862664687972017,  1.1721513422464978},
	{ 0.0,                 1.1721513422464978,  3.0862664687972017},
	{ 2.22474487139,      -2.22474487139,      -1.0},
	{ 2.22474487139,      -2.22474487139,       1.0},
	{ 1.1721513422464978, -3.0862664687972017,  0.0},
	{ 3.0862664687972017, -1.1721513422464978,  0.0},
	{ 2.22474487139,      -1.0,                -2.22474487139},
	{ 2.22474487139,       1.0,                -2.22474487139},
	{ 3.0862664687972017,  0.0,                -1.1721513422464978},
	{ 1.1721513422464978,  0.0,                -3.0862664687972017},
	{ 2.22474487139,      -1.0,                 2.22474487139},
	{ 2.22474487139,       1.0,                 2.22474487139},
	{ 1.1721513422464978,  0.0,                 3.0862664687972017},
	{ 3.0862664687972017,  0.0,                 1.1721513422464978},
	{ 2.22474487139,       2.22474487139,      -1.0},
	{ 2.22474487139,       2.22474487139,       1.0},
	{ 3.0862664687972017,  1.1721513422464978,  0.0},
	{ 1.1721513422464978,  3.0862664687972017,  0.0}
};

static struct Grad4 grad4[] = {
	{-0.753341017856078,    -0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624},
	{-0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301,    0.12128480194602098},
	{-0.7821684431180708,   -0.4321472685365301,    0.12128480194602098,  -0.4321472685365301},
	{-0.7821684431180708,    0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301},
	{-0.8586508742123365,   -0.508629699630796,     0.044802370851755174,  0.044802370851755174},
	{-0.8586508742123365,    0.044802370851755174, -0.508629699630796,     0.044802370851755174},
	{-0.8586508742123365,    0.044802370851755174,  0.044802370851755174, -0.508629699630796},
	{-0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842},
	{-0.37968289875261624,  -0.753341017856078,    -0.37968289875261624,  -0.37968289875261624},
	{-0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301,    0.12128480194602098},
	{-0.4321472685365301,   -0.7821684431180708,    0.12128480194602098,  -0.4321472685365301},
	{ 0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301},
	{-0.508629699630796,    -0.8586508742123365,    0.044802370851755174,  0.044802370851755174},
	{ 0.044802370851755174, -0.8586508742123365,   -0.508629699630796,     0.044802370851755174},
	{ 0.044802370851755174, -0.8586508742123365,    0.044802370851755174, -0.508629699630796},
	{-0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842},
	{-0.37968289875261624,  -0.37968289875261624,  -0.753341017856078,    -0.37968289875261624},
	{-0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708,    0.12128480194602098},
	{-0.4321472685365301,    0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301},
	{ 0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301},
	{-0.508629699630796,     0.044802370851755174, -0.8586508742123365,    0.044802370851755174},
	{ 0.044802370851755174, -0.508629699630796,    -0.8586508742123365,    0.044802370851755174},
	{ 0.044802370851755174,  0.044802370851755174, -0.8586508742123365,   -0.508629699630796},
	{-0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842},
	{-0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624,  -0.753341017856078},
	{-0.4321472685365301,   -0.4321472685365301,    0.12128480194602098,  -0.7821684431180708},
	{-0.4321472685365301,    0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708},
	{ 0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708},
	{-0.508629699630796,     0.044802370851755174,  0.044802370851755174, -0.8586508742123365},
	{ 0.044802370851755174, -0.508629699630796,     0.044802370851755174, -0.8586508742123365},
	{ 0.044802370851755174,  0.044802370851755174, -0.508629699630796,    -0.8586508742123365},
	{-0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062},
	{-0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537,    0.5794684678643381},
	{-0.7504883828755602,   -0.4004672082940195,    0.15296486218853164,   0.5029860367700724},
	{-0.7504883828755602,    0.15296486218853164,  -0.4004672082940195,    0.5029860367700724},
	{-0.8828161875373585,    0.08164729285680945,   0.08164729285680945,   0.4553054119602712},
	{-0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945,   0.8828161875373585},
	{-0.5029860367700724,   -0.15296486218853164,   0.4004672082940195,    0.7504883828755602},
	{-0.5029860367700724,    0.4004672082940195,   -0.15296486218853164,   0.7504883828755602},
	{-0.5794684678643381,    0.3239847771997537,    0.3239847771997537,    0.6740059517812944},
	{-0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537,    0.5794684678643381},
	{-0.4004672082940195,   -0.7504883828755602,    0.15296486218853164,   0.5029860367700724},
	{ 0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195,    0.5029860367700724},
	{ 0.08164729285680945,  -0.8828161875373585,    0.08164729285680945,   0.4553054119602712},
	{-0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945,   0.8828161875373585},
	{-0.15296486218853164,  -0.5029860367700724,    0.4004672082940195,    0.7504883828755602},
	{ 0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164,   0.7504883828755602},
	{ 0.3239847771997537,   -0.5794684678643381,    0.3239847771997537,    0.6740059517812944},
	{-0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944,    0.5794684678643381},
	{-0.4004672082940195,    0.15296486218853164,  -0.7504883828755602,    0.5029860367700724},
	{ 0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602,    0.5029860367700724},
	{ 0.08164729285680945,   0.08164729285680945,  -0.8828161875373585,    0.4553054119602712},
	{-0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712,    0.8828161875373585},
	{-0.15296486218853164,   0.4004672082940195,   -0.5029860367700724,    0.7504883828755602},
	{ 0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724,    0.7504883828755602},
	{ 0.3239847771997537,    0.3239847771997537,   -0.5794684678643381,    0.6740059517812944},
	{-0.6740059517812944,   -0.3239847771997537,    0.5794684678643381,   -0.3239847771997537},
	{-0.7504883828755602,   -0.4004672082940195,    0.5029860367700724,    0.15296486218853164},
	{-0.7504883828755602,    0.15296486218853164,   0.5029860367700724,   -0.4004672082940195},
	{-0.8828161875373585,    0.08164729285680945,   0.4553054119602712,    0.08164729285680945},
	{-0.4553054119602712,   -0.08164729285680945,   0.8828161875373585,   -0.08164729285680945},
	{-0.5029860367700724,   -0.15296486218853164,   0.7504883828755602,    0.4004672082940195},
	{-0.5029860367700724,    0.4004672082940195,    0.7504883828755602,   -0.15296486218853164},
	{-0.5794684678643381,    0.3239847771997537,    0.6740059517812944,    0.3239847771997537},
	{-0.3239847771997537,   -0.6740059517812944,    0.5794684678643381,   -0.3239847771997537},
	{-0.4004672082940195,   -0.7504883828755602,    0.5029860367700724,    0.15296486218853164},
	{ 0.15296486218853164,  -0.7504883828755602,    0.5029860367700724,   -0.4004672082940195},
	{ 0.08164729285680945,  -0.8828161875373585,    0.4553054119602712,    0.08164729285680945},
	{-0.08164729285680945,  -0.4553054119602712,    0.8828161875373585,   -0.08164729285680945},
	{-0.15296486218853164,  -0.5029860367700724,    0.7504883828755602,    0.4004672082940195},
	{ 0.4004672082940195,   -0.5029860367700724,    0.7504883828755602,   -0.15296486218853164},
	{ 0.3239847771997537,   -0.5794684678643381,    0.6740059517812944,    0.3239847771997537},
	{-0.3239847771997537,   -0.3239847771997537,    0.5794684678643381,   -0.6740059517812944},
	{-0.4004672082940195,    0.15296486218853164,   0.5029860367700724,   -0.7504883828755602},
	{ 0.15296486218853164,  -0.4004672082940195,    0.5029860367700724,   -0.7504883828755602},
	{ 0.08164729285680945,   0.08164729285680945,   0.4553054119602712,   -0.8828161875373585},
	{-0.08164729285680945,  -0.08164729285680945,   0.8828161875373585,   -0.4553054119602712},
	{-0.15296486218853164,   0.4004672082940195,    0.7504883828755602,   -0.5029860367700724},
	{ 0.4004672082940195,   -0.15296486218853164,   0.7504883828755602,   -0.5029860367700724},
	{ 0.3239847771997537,    0.3239847771997537,    0.6740059517812944,   -0.5794684678643381},
	{-0.6740059517812944,    0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537},
	{-0.7504883828755602,    0.5029860367700724,   -0.4004672082940195,    0.15296486218853164},
	{-0.7504883828755602,    0.5029860367700724,    0.15296486218853164,  -0.4004672082940195},
	{-0.8828161875373585,    0.4553054119602712,    0.08164729285680945,   0.08164729285680945},
	{-0.4553054119602712,    0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945},
	{-0.5029860367700724,    0.7504883828755602,   -0.15296486218853164,   0.4004672082940195},
	{-0.5029860367700724,    0.7504883828755602,    0.4004672082940195,   -0.15296486218853164},
	{-0.5794684678643381,    0.6740059517812944,    0.3239847771997537,    0.3239847771997537},
	{-0.3239847771997537,    0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537},
	{-0.4004672082940195,    0.5029860367700724,   -0.7504883828755602,    0.15296486218853164},
	{ 0.15296486218853164,   0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195},
	{ 0.08164729285680945,   0.4553054119602712,   -0.8828161875373585,    0.08164729285680945},
	{-0.08164729285680945,   0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945},
	{-0.15296486218853164,   0.7504883828755602,   -0.5029860367700724,    0.4004672082940195},
	{ 0.4004672082940195,    0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164},
	{ 0.3239847771997537,    0.6740059517812944,   -0.5794684678643381,    0.3239847771997537},
	{-0.3239847771997537,    0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944},
	{-0.4004672082940195,    0.5029860367700724,    0.15296486218853164,  -0.7504883828755602},
	{ 0.15296486218853164,   0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602},
	{ 0.08164729285680945,   0.4553054119602712,    0.08164729285680945,  -0.8828161875373585},
	{-0.08164729285680945,   0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712},
	{-0.15296486218853164,   0.7504883828755602,    0.4004672082940195,   -0.5029860367700724},
	{ 0.4004672082940195,    0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724},
	{ 0.3239847771997537,    0.6740059517812944,    0.3239847771997537,   -0.5794684678643381},
	{ 0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537},
	{ 0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195,    0.15296486218853164},
	{ 0.5029860367700724,   -0.7504883828755602,    0.15296486218853164,  -0.4004672082940195},
	{ 0.4553054119602712,   -0.8828161875373585,    0.08164729285680945,   0.08164729285680945},
	{ 0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945},
	{ 0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164,   0.4004672082940195},
	{ 0.7504883828755602,   -0.5029860367700724,    0.4004672082940195,   -0.15296486218853164},
	{ 0.6740059517812944,   -0.5794684678643381,    0.3239847771997537,    0.3239847771997537},
	{ 0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537},
	{ 0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602,    0.15296486218853164},
	{ 0.5029860367700724,    0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195},
	{ 0.4553054119602712,    0.08164729285680945,  -0.8828161875373585,    0.08164729285680945},
	{ 0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945},
	{ 0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724,    0.4004672082940195},
	{ 0.7504883828755602,    0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164},
	{ 0.6740059517812944,    0.3239847771997537,   -0.5794684678643381,    0.3239847771997537},
	{ 0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944},
	{ 0.5029860367700724,   -0.4004672082940195,    0.15296486218853164,  -0.7504883828755602},
	{ 0.5029860367700724,    0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602},
	{ 0.4553054119602712,    0.08164729285680945,   0.08164729285680945,  -0.8828161875373585},
	{ 0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712},
	{ 0.7504883828755602,   -0.15296486218853164,   0.4004672082940195,   -0.5029860367700724},
	{ 0.7504883828755602,    0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724},
	{ 0.6740059517812944,    0.3239847771997537,    0.3239847771997537,   -0.5794684678643381},
	{ 0.03381941603233842,   0.03381941603233842,   0.03381941603233842,   0.9982828964265062},
	{-0.044802370851755174, -0.044802370851755174,  0.508629699630796,     0.8586508742123365},
	{-0.044802370851755174,  0.508629699630796,    -0.044802370851755174,  0.8586508742123365},
	{-0.12128480194602098,   0.4321472685365301,    0.4321472685365301,    0.7821684431180708},
	{ 0.508629699630796,    -0.044802370851755174, -0.044802370851755174,  0.8586508742123365},
	{ 0.4321472685365301,   -0.12128480194602098,   0.4321472685365301,    0.7821684431180708},
	{ 0.4321472685365301,    0.4321472685365301,   -0.12128480194602098,   0.7821684431180708},
	{ 0.37968289875261624,   0.37968289875261624,   0.37968289875261624,   0.753341017856078},
	{ 0.03381941603233842,   0.03381941603233842,   0.9982828964265062,    0.03381941603233842},
	{-0.044802370851755174,  0.044802370851755174,  0.8586508742123365,    0.508629699630796},
	{-0.044802370851755174,  0.508629699630796,     0.8586508742123365,   -0.044802370851755174},
	{-0.12128480194602098,   0.4321472685365301,    0.7821684431180708,    0.4321472685365301},
	{ 0.508629699630796,    -0.044802370851755174,  0.8586508742123365,   -0.044802370851755174},
	{ 0.4321472685365301,   -0.12128480194602098,   0.7821684431180708,    0.4321472685365301},
	{ 0.4321472685365301,    0.4321472685365301,    0.7821684431180708,   -0.12128480194602098},
	{ 0.37968289875261624,   0.37968289875261624,   0.753341017856078,     0.37968289875261624},
	{ 0.03381941603233842,   0.9982828964265062,    0.03381941603233842,   0.03381941603233842},
	{-0.044802370851755174,  0.8586508742123365,   -0.044802370851755174,  0.508629699630796},
	{-0.044802370851755174,  0.8586508742123365,    0.508629699630796,    -0.044802370851755174},
	{-0.12128480194602098,   0.7821684431180708,    0.4321472685365301,    0.4321472685365301},
	{ 0.508629699630796,     0.8586508742123365,   -0.044802370851755174, -0.044802370851755174},
	{ 0.4321472685365301,    0.7821684431180708,   -0.12128480194602098,   0.4321472685365301},
	{ 0.4321472685365301,    0.7821684431180708,    0.4321472685365301,   -0.12128480194602098},
	{ 0.37968289875261624,   0.753341017856078,     0.37968289875261624,   0.37968289875261624},
	{ 0.9982828964265062,    0.03381941603233842,   0.03381941603233842,   0.03381941603233842},
	{ 0.8586508742123365,   -0.044802370851755174, -0.044802370851755174,  0.508629699630796},
	{ 0.8586508742123365,   -0.044802370851755174,  0.508629699630796,    -0.044802370851755174},
	{ 0.7821684431180708,   -0.12128480194602098,   0.4321472685365301,    0.4321472685365301},
	{ 0.8586508742123365,    0.508629699630796,    -0.044802370851755174, -0.044802370851755174},
	{ 0.7821684431180708,    0.4321472685365301,   -0.12128480194602098,   0.4321472685365301},
	{ 0.7821684431180708,    0.4321472685365301,    0.4321472685365301,   -0.12128480194602098},
	{ 0.753341017856078,     0.37968289875261624,   0.37968289875261624,   0.37968289875261624}
};


static void setup_gradients(void)
{
	static int already_did = 0;
	int i;

	if (already_did)
		return;
	already_did = 1;

	for (i = 0; (size_t) i < ARRAYSIZE(grad2); i++) {
		grad2[i].dx /= N2; grad2[i].dy /= N2;
	}
	for (int i = 0; i < PSIZE; i++) {
		GRADIENTS_2D[i] = grad2[i % ARRAYSIZE(grad2)];
	}

	for (i = 0; (size_t) i < ARRAYSIZE(grad3); i++) {
		grad3[i].dx /= N3; grad3[i].dy /= N3; grad3[i].dz /= N3;
	}
	for (i = 0; i < PSIZE; i++) {
		GRADIENTS_3D[i] = grad3[i % ARRAYSIZE(grad3)];
	}

	for (i = 0; (size_t) i < ARRAYSIZE(grad4); i++) {
		grad4[i].dx /= N4; grad4[i].dy /= N4; grad4[i].dz /= N4; grad4[i].dw /= N4;
	}
	for (i = 0; i < PSIZE; i++) {
		GRADIENTS_4D[i] = grad4[i % ARRAYSIZE(grad4)];
	}
}

static void setup_lattice_points(void)
{
	static int already_did = 0;
	int i;

	if (already_did)
		return;
	already_did = 1;

	LOOKUP_2D[0] = new_LatticePoint2D(1, 0);
	LOOKUP_2D[1] = new_LatticePoint2D(0, 0);
	LOOKUP_2D[2] = new_LatticePoint2D(1, 1);
	LOOKUP_2D[3] = new_LatticePoint2D(0, 1);

	for (i = 0; i < 8; i++) {
		int i1, j1, k1, i2, j2, k2;
		i1 = (i >> 0) & 1; j1 = (i >> 1) & 1; k1 = (i >> 2) & 1;
		i2 = i1 ^ 1; j2 = j1 ^ 1; k2 = k1 ^ 1;

		// The two points within this octant, one from each of the two cubic half-lattices.
		struct LatticePoint3D *c0 = new_LatticePoint3D(i1, j1, k1, 0);
		struct LatticePoint3D *c1 = new_LatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1);

		// Each single step away on the first half-lattice.
		struct LatticePoint3D *c2 = new_LatticePoint3D(i1 ^ 1, j1, k1, 0);
		struct LatticePoint3D *c3 = new_LatticePoint3D(i1, j1 ^ 1, k1, 0);
		struct LatticePoint3D *c4 = new_LatticePoint3D(i1, j1, k1 ^ 1, 0);

		// Each single step away on the second half-lattice.
		struct LatticePoint3D *c5 = new_LatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1);
		struct LatticePoint3D *c6 = new_LatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1);
		struct LatticePoint3D *c7 = new_LatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1);

		// First two are guaranteed.
		c0->nextOnFailure = c0->nextOnSuccess = c1;
		c1->nextOnFailure = c1->nextOnSuccess = c2;

		// Once we find one on the first half-lattice, the rest are out.
		// In addition, knowing c2 rules out c5.
		c2->nextOnFailure = c3; c2->nextOnSuccess = c6;
		c3->nextOnFailure = c4; c3->nextOnSuccess = c5;
		c4->nextOnFailure = c4->nextOnSuccess = c5;

		// Once we find one on the second half-lattice, the rest are out.
		c5->nextOnFailure = c6; c5->nextOnSuccess = NULL;
		c6->nextOnFailure = c7; c6->nextOnSuccess = NULL;
		c7->nextOnFailure = c7->nextOnSuccess = NULL;

		LOOKUP_3D[i] = c0;
	}

	for (i = 0; i < 16; i++) {
		VERTICES_4D[i] = new_LatticePoint4D((i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1, (i >> 3) & 1);
	}
}

/* Free up all the LatticePoint stuff allocated in setup_lattice_points() */
void OpenSimplex2F_shutdown(void)
{
	int i;
	struct LatticePoint3D *list[8];
	int count = 0;

	for (i = 0; (size_t) i < ARRAYSIZE(LOOKUP_2D); i++) {
		if (LOOKUP_2D[i]) {
			free(LOOKUP_2D[i]);
			LOOKUP_2D[i] = NULL;
		}
	}

	for (i = 0; (size_t) i < ARRAYSIZE(LOOKUP_3D); i++) {
		if (LOOKUP_3D[i]) {
			/* There are cycles in the tree LOOKUP_3D[i], so we have to
			 * make a list of unique pointers within the tree to free them. 
			 * We know there are only 8 unique pointers in the tree.
			 */
			count = 0;
			memset(list, 0, sizeof(list));
			find_unique_pointers(LOOKUP_3D[i], list, &count);
			free_LatticePoint3D(list, count);
			LOOKUP_3D[i] = NULL;
		}
	}

	for (i = 0; (size_t) i < ARRAYSIZE(VERTICES_4D); i++) {
		if (VERTICES_4D[i]) {
			free(VERTICES_4D[i]);
			VERTICES_4D[i] = NULL;
		}
	}
}

void OpenSimplex2F_free(struct OpenSimplex2F_context *ctx)
{
	if (ctx->perm)
		free(ctx->perm);
	if (ctx->permGrad2)
		free(ctx->permGrad2);
	if (ctx->permGrad3)
		free(ctx->permGrad3);
	if (ctx->permGrad4)
		free(ctx->permGrad4);
	free(ctx);
}

int OpenSimplex2F(int64_t seed, struct OpenSimplex2F_context **ctx)
{
	struct OpenSimplex2F_context *c;
	int i, *source;

	setup_gradients();
	setup_lattice_points();

	source = calloc(1, sizeof(*source) * PSIZE);
	if (!source)
		return -1;

	c = calloc(1, sizeof(**ctx));
	if (!c) {
		free(source);
		return -1;
	}
	*ctx = c;

	c->perm = calloc(1, sizeof(*c->perm) * PSIZE);
	c->permGrad2 = calloc(1, sizeof(*c->permGrad2) * PSIZE);
	c->permGrad3 = calloc(1, sizeof(*c->permGrad3) * PSIZE);
	c->permGrad4 = calloc(1, sizeof(*c->permGrad4) * PSIZE);

	if (!c->perm || !c->permGrad2 || !c->permGrad3 || !c->permGrad4) {
		OpenSimplex2F_free(*ctx);
		*ctx = NULL;
		free(source);
		return -1;
	}

	for (i = 0; i < PSIZE; i++)
		source[i] = (int16_t) i;

	for (int i = PSIZE - 1; i >= 0; i--) {
		seed = (uint64_t) seed * (uint64_t) 6364136223846793005LL + (uint64_t) 1442695040888963407LL;
		int r = (int)((seed + 31) % (i + 1));
		if (r < 0)
			r += (i + 1);
		c->perm[i] = source[r];
		c->permGrad2[i] = GRADIENTS_2D[c->perm[i]];
		c->permGrad3[i] = GRADIENTS_3D[c->perm[i]];
		c->permGrad4[i] = GRADIENTS_4D[c->perm[i]];
		source[r] = source[i];
	}
	free(source);
	return 0;
}

static int fastFloor(double x)
{
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
}

/*
 * 2D Simplex noise base.
 * Lookup table implementation inspired by DigitalShadow.
 */
static double noise2_Base(struct OpenSimplex2F_context *ctx, double xs, double ys)
{
	double value = 0;
	int i;

	// Get base points and offsets
	int xsb = fastFloor(xs), ysb = fastFloor(ys);
	double xsi = xs - xsb, ysi = ys - ysb;

	// Index to point list
	int index = (int)((ysi - xsi) / 2 + 1);

	double ssi = (xsi + ysi) * -0.211324865405187;
	double xi = xsi + ssi, yi = ysi + ssi;

	// Point contributions
	for (i = 0; i < 3; i++) {
		struct LatticePoint2D c = *(LOOKUP_2D[index + i]);

		double dx = xi + c.dx, dy = yi + c.dy;
		double attn = 0.5 - dx * dx - dy * dy;
		if (attn <= 0) continue;

		int pxm = (xsb + c.xsv) & PMASK, pym = (ysb + c.ysv) & PMASK;
		struct Grad2 grad = ctx->permGrad2[ctx->perm[pxm] ^ pym];
		double extrapolation = grad.dx * dx + grad.dy * dy;

		attn *= attn;
		value += attn * attn * extrapolation;
	}

	return value;
}

/**
 * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
 * Lookup table implementation inspired by DigitalShadow.
 * It was actually faster to narrow down the points in the loop itself,
 * than to build up the index with enough info to isolate 4 points.
 */
static double noise3_BCC(struct OpenSimplex2F_context *ctx, double xr, double yr, double zr)
{
	// Get base and offsets inside cube of first lattice.
	int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
	double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

	// Identify which octant of the cube we're in. This determines which cell
	// in the other cubic lattice we're in, and also narrows down one point on each.
	int xht = (int)(xri + 0.5), yht = (int)(yri + 0.5), zht = (int)(zri + 0.5);
	int index = (xht << 0) | (yht << 1) | (zht << 2);

	// Point contributions
	double value = 0;
	struct LatticePoint3D *c = LOOKUP_3D[index];
	while (c != NULL) {
		double dxr = xri + c->dxr, dyr = yri + c->dyr, dzr = zri + c->dzr;
		double attn = 0.5 - dxr * dxr - dyr * dyr - dzr * dzr;
		if (attn < 0) {
			c = c->nextOnFailure;
		} else {
			int pxm = (xrb + c->xrv) & PMASK, pym = (yrb + c->yrv) & PMASK, pzm = (zrb + c->zrv) & PMASK;
			struct Grad3 grad = ctx->permGrad3[ctx->perm[ctx->perm[pxm] ^ pym] ^ pzm];
			double extrapolation = grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

			attn *= attn;
			value += attn * attn * extrapolation;
			c = c->nextOnSuccess;
		}
	}
	return value;
}

/**
 * 2D Simplex noise, standard lattice orientation.
 */
double OpenSimplex2F_noise2(struct OpenSimplex2F_context *ctx, double x, double y)
{
	// Get points for A2* lattice
	double s = 0.366025403784439 * (x + y);
	double xs = x + s, ys = y + s;

	return noise2_Base(ctx, xs, ys);
}

/**
 * 2D Simplex noise, with Y pointing down the main diagonal.
 * Might be better for a 2D sandbox style game, where Y is vertical.
 * Probably slightly less optimal for heightmaps or continent maps.
 */
double OpenSimplex2F_noise2_XBeforeY(struct OpenSimplex2F_context *ctx, double x, double y)
{
	// Skew transform and rotation baked into one.
	double xx = x * 0.7071067811865476;
	double yy = y * 1.224744871380249;

	return noise2_Base(ctx, yy + xx, yy - xx);
}

/**
 * 3D Re-oriented 4-point BCC noise, classic orientation.
 * Proper substitute for 3D Simplex in light of Forbidden Formulae.
 * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
 */
double OpenSimplex2F_noise3_Classic(struct OpenSimplex2F_context *ctx, double x, double y, double z)
{
	// Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
	// If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
	// Orthonormal rotation. Not a skew transform.
	double r = (2.0 / 3.0) * (x + y + z);
	double xr = r - x, yr = r - y, zr = r - z;

	// Evaluate both lattices to form a BCC lattice.
	return noise3_BCC(ctx, xr, yr, zr);
}

/**
 * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Y).
 * Recommended for 3D terrain and time-varied animations.
 * The Z coordinate should always be the "different" coordinate in your use case.
 * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
 * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
 * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
 */
double OpenSimplex2F_noise3_XYBeforeZ(struct OpenSimplex2F_context *ctx, double x, double y, double z)
{
	// Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
	// Orthonormal rotation. Not a skew transform.
	double xy = x + y;
	double s2 = xy * -0.211324865405187;
	double zz = z * 0.577350269189626;
	double xr = x + s2 - zz, yr = y + s2 - zz;
	double zr = xy * 0.577350269189626 + zz;

	// Evaluate both lattices to form a BCC lattice.
	return noise3_BCC(ctx, xr, yr, zr);
}

/**
 * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Z).
 * Recommended for 3D terrain and time-varied animations.
 * The Y coordinate should always be the "different" coordinate in your use case.
 * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
 * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
 * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
 */
double OpenSimplex2F_noise3_XZBeforeY(struct OpenSimplex2F_context *ctx, double x, double y, double z)
{
	// Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
	// Orthonormal rotation. Not a skew transform.
	double xz = x + z;
	double s2 = xz * -0.211324865405187;
	double yy = y * 0.577350269189626;
	double xr = x + s2 - yy; double zr = z + s2 - yy;
	double yr = xz * 0.577350269189626 + yy;

	// Evaluate both lattices to form a BCC lattice.
	return noise3_BCC(ctx, xr, yr, zr);
}

/**
 * 4D OpenSimplex2F noise base.
 * Current implementation not fully optimized by lookup tables.
 * But still comes out slightly ahead of Gustavson's Simplex in tests.
 */
static double noise4_Base(struct OpenSimplex2F_context *ctx, double xs, double ys, double zs, double ws)
{
	double value = 0;

	// Get base points and offsets
	int xsb = fastFloor(xs), ysb = fastFloor(ys), zsb = fastFloor(zs), wsb = fastFloor(ws);
	double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

	// If we're in the lower half, flip so we can repeat the code for the upper half. We'll flip back later.
	double siSum = xsi + ysi + zsi + wsi;
	double ssi = siSum * 0.309016994374947; // Prep for vertex contributions.
	bool inLowerHalf = (siSum < 2);
	if (inLowerHalf) {
		xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
		siSum = 4 - siSum;
	}

	// Consider opposing vertex pairs of the octahedron formed by the central cross-section of the stretched tesseract
	double aabb = xsi + ysi - zsi - wsi, abab = xsi - ysi + zsi - wsi, abba = xsi - ysi - zsi + wsi;
	double aabbScore = fabs(aabb), ababScore = fabs(abab), abbaScore = fabs(abba);

	// Find the closest point on the stretched tesseract as if it were the upper half
	int vertexIndex, via, vib;
	double asi, bsi;
	if (aabbScore > ababScore && aabbScore > abbaScore) {
		if (aabb > 0) {
			asi = zsi; bsi = wsi; vertexIndex = B0011; via = B0111; vib = B1011;
		} else {
			asi = xsi; bsi = ysi; vertexIndex = B1100; via = B1101; vib = B1110;
		}
	} else if (ababScore > abbaScore) {
		if (abab > 0) {
			asi = ysi; bsi = wsi; vertexIndex = B0101; via = B0111; vib = B1101;
		} else {
			asi = xsi; bsi = zsi; vertexIndex = B1010; via = B1011; vib = B1110;
		}
	} else {
		if (abba > 0) {
			asi = ysi; bsi = zsi; vertexIndex = B1001; via = B1011; vib = B1101;
		} else {
			asi = xsi; bsi = wsi; vertexIndex = B0110; via = B0111; vib = B1110;
		}
	}
	if (bsi > asi) {
		via = vib;
		double temp = bsi;
		bsi = asi;
		asi = temp;
	}
	if (siSum + asi > 3) {
		vertexIndex = via;
		if (siSum + bsi > 4) {
			vertexIndex = B1111;
		}
	}

	// Now flip back if we're actually in the lower half.
	if (inLowerHalf) {
		xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
		vertexIndex ^= B1111;
	}

	// Five points to add, total, from five copies of the A4 lattice.
	for (int i = 0; i < 5; i++) {

		// Update xsb/etc. and add the lattice point's contribution.
		struct LatticePoint4D c = *(VERTICES_4D[vertexIndex]);
		xsb += c.xsv; ysb += c.ysv; zsb += c.zsv; wsb += c.wsv;
		double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;
		double dx = xi + c.dx, dy = yi + c.dy, dz = zi + c.dz, dw = wi + c.dw;
		double attn = 0.5 - dx * dx - dy * dy - dz * dz - dw * dw;
		if (attn > 0) {
			int pxm = xsb & PMASK, pym = ysb & PMASK, pzm = zsb & PMASK, pwm = wsb & PMASK;
			struct Grad4 grad = ctx->permGrad4[ctx->perm[ctx->perm[ctx->perm[pxm] ^ pym] ^ pzm] ^ pwm];
			double ramped = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

			attn *= attn;
			value += attn * attn * ramped;
		}

		// Maybe this helps the compiler/JVM/LLVM/etc. know we can end the loop here. Maybe not.
		if (i == 4) break;

		// Update the relative skewed coordinates to reference the vertex we just added.
		// Rather, reference its counterpart on the lattice copy that is shifted down by
		// the vector <-0.2, -0.2, -0.2, -0.2>
		xsi += c.xsi; ysi += c.ysi; zsi += c.zsi; wsi += c.wsi;
		ssi += c.ssiDelta;

		// Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
		double score0 = 1.0 + ssi * (-1.0 / 0.309016994374947); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
		vertexIndex = B0000;
		if (xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0) {
			vertexIndex = B0001;
		}
		else if (ysi > xsi && ysi >= zsi && ysi >= wsi && ysi >= score0) {
			vertexIndex = B0010;
		}
		else if (zsi > xsi && zsi > ysi && zsi >= wsi && zsi >= score0) {
			vertexIndex = B0100;
		}
		else if (wsi > xsi && wsi > ysi && wsi > zsi && wsi >= score0) {
			vertexIndex = B1000;
		}
	}

	return value;
}

/**
 * 4D OpenSimplex2F noise, classic lattice orientation.
 */
 double OpenSimplex2F_noise4_Classic(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w)
{
	// Get points for A4 lattice
	double s = -0.138196601125011 * (x + y + z + w);
	double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

	return noise4_Base(ctx, xs, ys, zs, ws);
}

/**
 * 4D OpenSimplex2F noise, with XY and ZW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
 * Recommended for noise(x, y, sin(time), cos(time)) trick.
 */
double OpenSimplex2F_noise4_XYBeforeZW(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w)
{

	double s2 = (x + y) * -0.178275657951399372 + (z + w) * 0.215623393288842828;
	double t2 = (z + w) * -0.403949762580207112 + (x + y) * -0.375199083010075342;
	double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

	return noise4_Base(ctx, xs, ys, zs, ws);
}

/**
 * 4D OpenSimplex2F noise, with XZ and YW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
 */
double OpenSimplex2F_noise4_XZBeforeYW(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w)
{
	double s2 = (x + z) * -0.178275657951399372 + (y + w) * 0.215623393288842828;
	double t2 = (y + w) * -0.403949762580207112 + (x + z) * -0.375199083010075342;
	double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

	return noise4_Base(ctx, xs, ys, zs, ws);
}

/**
 * 4D OpenSimplex2F noise, with XYZ oriented like noise3_Classic,
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 */
double OpenSimplex2F_noise4_XYZBeforeW(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w)
{
	double xyz = x + y + z;
	double ww = w * 0.2236067977499788;
	double s2 = xyz * -0.16666666666666666 + ww;
	double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

	return noise4_Base(ctx, xs, ys, zs, ws);
}

