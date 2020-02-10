/*
 * Noise Demo 2D.
 * Replace OpenSimplex2S with OpenSimplex2F to demo OpenSimplex2F instead.
 */

import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import java.io.*;
import javax.swing.*;

public class NoiseDemoEvalOnly
{
	private static final int WIDTH = 1024;
	private static final int HEIGHT = 1024;
	private static final double PERIOD = 64.0;
	private static final int OFF_X = 2048;
	private static final int OFF_Y = 2048;
	
	private static final double FREQ = 1.0 / PERIOD;

	public static void main(String[] args)
			throws IOException {
		
		// Initialize
		OpenSimplex2S noise = new OpenSimplex2S(1);
		
		// Image
		BufferedImage image = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
		for (int y = 0; y < HEIGHT; y++)
		{
			for (int x = 0; x < WIDTH; x++)
			{
				double value = noise.noise3_XZBeforeY((x + OFF_X) * FREQ, 0.0, (y + OFF_Y) * FREQ);
				
				//if (value < -1) value = -1;
				//else if (value > 1) value = 1;
				
				int rgb = 0x010101 * (int)((value + 1) * 127.5);
				image.setRGB(x, y, rgb);
			}
		}
		
		// Save it or show it
		if (args.length > 0 && args[0] != null) {
			ImageIO.write(image, "png", new File(args[0]));
			System.out.println("Saved image as " + args[0]);
		} else {
			JFrame frame = new JFrame();
			JLabel imageLabel = new JLabel();
			imageLabel.setIcon(new ImageIcon(image));
			frame.add(imageLabel);
			frame.pack();
			frame.setResizable(false);
			frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
			frame.setVisible(true);
		}
		
	}
}