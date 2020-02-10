import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

class NoiseMetrics2_FSS {
	
	static final int N_PREP_ITERATIONS = 32;
	static final int N_TIMED_ITERATIONS = 256;
	
	static final int WIDTH = 2048;
	static final int HEIGHT = 2048;
	static final double NOISE_EVAL_PERIOD = 64.0;
	
	static final double NOISE_EVAL_FREQ = 1.0 / NOISE_EVAL_PERIOD;
	
	public static void main(String[] args) {
		
		List<NoiseTimer> noiseTimers = new ArrayList<>();
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "OpenSimplex2F Area Generation"; }
			OpenSimplex2F.GenerateContext2D ctx = new OpenSimplex2F.GenerateContext2D(NOISE_EVAL_FREQ, NOISE_EVAL_FREQ, 1.0);
			OpenSimplex2F noise = new OpenSimplex2F(0);
			
			void test(int offX, int offY) {
				double[][] buffer = new double[HEIGHT][WIDTH];
				noise.generate2(ctx, buffer, offX, offY, WIDTH, HEIGHT, 0, 0);
			}
			
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "OpenSimplex2F Evaluation"; }
			OpenSimplex2F noise = new OpenSimplex2F(0);
			
			void test(int offX, int offY) {
				double[][] buffer = new double[HEIGHT][WIDTH];
				for (int y = 0; y < HEIGHT; y++) {
					for (int x = 0; x < WIDTH; x++) {
						buffer[y][x] = noise.noise2((x + offX) * NOISE_EVAL_FREQ, (y + offY) * NOISE_EVAL_FREQ);
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "Gustavson's Simplex Noise Implementation"; }
			
			void test(int offX, int offY) {
				double[][] buffer = new double[HEIGHT][WIDTH];
				for (int y = 0; y < HEIGHT; y++) {
					for (int x = 0; x < WIDTH; x++) {
						buffer[y][x] = SimplexNoise.noise((x + offX) * NOISE_EVAL_FREQ, (y + offY) * NOISE_EVAL_FREQ);
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "FastNoise Simplex Implementation"; }
			FastNoise fastNoise = new FastNoise(0);
			{ fastNoise.SetNoiseType(FastNoise.NoiseType.Simplex); }
			
			void test(int offX, int offY) {
				double[][] buffer = new double[HEIGHT][WIDTH];
				for (int y = 0; y < HEIGHT; y++) {
					for (int x = 0; x < WIDTH; x++) {
						buffer[y][x] = fastNoise.GetNoise((float)((x + offX) * NOISE_EVAL_FREQ), (float)((y + offY) * NOISE_EVAL_FREQ));
					}
				}
			}
		});
		
		noiseTimers.add(new NoiseTimer() {
			{ name = "FastNoise Oldperlin Implementation"; }
			FastNoise fastNoise = new FastNoise(0);
			{ fastNoise.SetNoiseType(FastNoise.NoiseType.Perlin); }
			
			void test(int offX, int offY) {
				double[][] buffer = new double[HEIGHT][WIDTH];
				for (int y = 0; y < HEIGHT; y++) {
					for (int x = 0; x < WIDTH; x++) {
						buffer[y][x] = fastNoise.GetNoise((float)((x + offX) * NOISE_EVAL_FREQ), (float)((y + offY) * NOISE_EVAL_FREQ));
					}
				}
			}
		});
		
		System.out.println("Number of prep iterations: " + N_PREP_ITERATIONS);
		System.out.println("Number of timed iterations: " + N_TIMED_ITERATIONS);
		System.out.println("Size: " + WIDTH  + "x" + HEIGHT);
		System.out.println("Noise Period: " + NOISE_EVAL_PERIOD);
		
		for (NoiseTimer timer : noiseTimers) {
			System.out.println();
			System.out.println("---- " + timer.name + " (No Image Display) ----");
			timer.run();
			System.out.println("Total milliseconds: " + timer.time);
			System.out.println("Nanoseconds per generated value: " + (timer.time * 1_000_000.0 / (N_TIMED_ITERATIONS * WIDTH * HEIGHT)));
		}
	
	}

	static abstract class NoiseTimer {
		String name;
		int time, sum;
		
		abstract void test(int offX, int offY);
		
		void run() {
			int counter = 0;
			for (int ie = 0; ie < N_PREP_ITERATIONS + N_TIMED_ITERATIONS; ie++) {
				long start = System.currentTimeMillis();
				test(counter * WIDTH, 0);
				long elapsed = System.currentTimeMillis() - start;
				
				//sum += Arrays.stream(values).sum();
				if (ie >= N_PREP_ITERATIONS) time += elapsed;
				
				counter++;
			}
		}
	}
	
}