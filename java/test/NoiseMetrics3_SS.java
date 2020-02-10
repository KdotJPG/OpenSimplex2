import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

class NoiseMetrics3_SS {
	
	static final int N_PREP_ITERATIONS = 32;
	static final int N_TIMED_ITERATIONS = 256;
	
	static final int WIDTH = 128;
	static final int HEIGHT = 128;
	static final int DEPTH = 128;
	static final double NOISE_EVAL_PERIOD = 64.0;
	
	static final double NOISE_EVAL_FREQ = 1.0 / NOISE_EVAL_PERIOD;
	
	public static void main(String[] args) {
		
		List<NoiseTimer> noiseTimers = new ArrayList<>();
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "OpenSimplex2S Area Generation Classic Orientation"; }
			OpenSimplex2S.GenerateContext3D ctx = new OpenSimplex2S.GenerateContext3D(OpenSimplex2S.LatticeOrientation3D.Classic, NOISE_EVAL_FREQ, NOISE_EVAL_FREQ, NOISE_EVAL_FREQ, 1.0);
			OpenSimplex2S noise = new OpenSimplex2S(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				noise.generate3(ctx, buffer, offX, offY, offZ, WIDTH, HEIGHT, DEPTH, 0, 0, 0);
			}
			
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "OpenSimplex2S Area Generation PlaneFirst Orientation"; }
			OpenSimplex2S.GenerateContext3D ctx = new OpenSimplex2S.GenerateContext3D(OpenSimplex2S.LatticeOrientation3D.PlaneFirst, NOISE_EVAL_FREQ, NOISE_EVAL_FREQ, NOISE_EVAL_FREQ, 1.0);
			OpenSimplex2S noise = new OpenSimplex2S(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				noise.generate3(ctx, buffer, offX, offY, offZ, WIDTH, HEIGHT, DEPTH, 0, 0, 0);
			}
			
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "OpenSimplex2S Evaluation Classic Orientation"; }
			OpenSimplex2S noise = new OpenSimplex2S(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				for (int z = 0; z < DEPTH; z++) {
					for (int y = 0; y < HEIGHT; y++) {
						for (int x = 0; x < WIDTH; x++) {
							buffer[z][y][x] = noise.noise3_Classic((x + offX) * NOISE_EVAL_FREQ, (y + offY) * NOISE_EVAL_FREQ, (z + offZ) * NOISE_EVAL_FREQ);
						}
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "OpenSimplex2S Evaluation PlaneFirst Orientation"; }
			OpenSimplex2S noise = new OpenSimplex2S(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				for (int z = 0; z < DEPTH; z++) {
					for (int y = 0; y < HEIGHT; y++) {
						for (int x = 0; x < WIDTH; x++) {
							buffer[z][y][x] = noise.noise3_PlaneFirst((x + offX) * NOISE_EVAL_FREQ, (y + offY) * NOISE_EVAL_FREQ, (z + offZ) * NOISE_EVAL_FREQ);
						}
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "DigitalShadow's Optimized OpenSimplex Noise Classic Orientation"; }
			OpenSimplex noise = new OpenSimplex(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				for (int z = 0; z < DEPTH; z++) {
					for (int y = 0; y < HEIGHT; y++) {
						for (int x = 0; x < WIDTH; x++) {
							buffer[z][y][x] = noise.eval((x + offX) * NOISE_EVAL_FREQ * 2, (y + offY) * NOISE_EVAL_FREQ * 2, (z + offZ) * NOISE_EVAL_FREQ * 2);
						}
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "DigitalShadow's Optimized OpenSimplex Noise XZBeforeY Orientation"; }
			OpenSimplex noise = new OpenSimplex(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				for (int z = 0; z < DEPTH; z++) {
					for (int y = 0; y < HEIGHT; y++) {
						for (int x = 0; x < WIDTH; x++) {
							buffer[z][y][x] = noise.eval3_XZBeforeY((x + offX) * NOISE_EVAL_FREQ * 2, (y + offY) * NOISE_EVAL_FREQ * 2, (z + offZ) * NOISE_EVAL_FREQ * 2);
						}
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "Unoptimized OpenSimplex Noise Classic Orientation"; }
			OpenSimplexUnoptimized noise = new OpenSimplexUnoptimized(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				for (int z = 0; z < DEPTH; z++) {
					for (int y = 0; y < HEIGHT; y++) {
						for (int x = 0; x < WIDTH; x++) {
							buffer[z][y][x] = noise.eval((x + offX) * NOISE_EVAL_FREQ * 2, (y + offY) * NOISE_EVAL_FREQ * 2, (z + offZ) * NOISE_EVAL_FREQ * 2);
						}
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "Legacy OpenSimplex Noise XZBeforeY Orientation"; }
			OpenSimplexUnoptimized noise = new OpenSimplexUnoptimized(0);
			
			void test(int offX, int offY, int offZ) {
				double[][][] buffer = new double[DEPTH][HEIGHT][WIDTH];
				for (int z = 0; z < DEPTH; z++) {
					for (int y = 0; y < HEIGHT; y++) {
						for (int x = 0; x < WIDTH; x++) {
							buffer[z][y][x] = noise.eval3_XZBeforeY((x + offX) * NOISE_EVAL_FREQ * 2, (y + offY) * NOISE_EVAL_FREQ * 2, (z + offZ) * NOISE_EVAL_FREQ * 2);
						}
					}
				}
			}
		});
		
		System.out.println("Number of prep iterations: " + N_PREP_ITERATIONS);
		System.out.println("Number of timed iterations: " + N_TIMED_ITERATIONS);
		System.out.println("Size: " + WIDTH  + "x" + HEIGHT + "x" + DEPTH);
		System.out.println("Noise Period: " + NOISE_EVAL_PERIOD);
		
		for (NoiseTimer timer : noiseTimers) {
			System.out.println();
			System.out.println("---- " + timer.name + " (No Image Display) ----");
			timer.run();
			System.out.println("Total milliseconds: " + timer.time);
			System.out.println("Nanoseconds per generated value: " + (timer.time * 1_000_000.0 / (N_TIMED_ITERATIONS * WIDTH * HEIGHT * DEPTH)));
		}
	
	}

	static abstract class NoiseTimer {
		String name;
		int time, sum;
		
		abstract void test(int offX, int offY, int offZ);
		
		void run() {
			int counter = 0;
			for (int ie = 0; ie < N_PREP_ITERATIONS + N_TIMED_ITERATIONS; ie++) {
				long start = System.currentTimeMillis();
				test(counter * WIDTH, 0, 0);
				long elapsed = System.currentTimeMillis() - start;
				
				//sum += Arrays.stream(values).sum();
				if (ie >= N_PREP_ITERATIONS) time += elapsed;
				
				counter++;
			}
		}
	}
	
}