

// BV Ue7 SS2023 Vorgabe
//
// Copyright (C) 2023 by Klaus Jung
// All rights reserved.
// Date: 2023-03-23 		   	  	  		
package bv_ss23;




public class DiscreteCosineTransform {
		
	
	
	  
    public void processDCT(RasterImage originalImage, RasterImage dctImage, RasterImage reconstructedImage, int numCoefficients) {
        int width = originalImage.width;
        int height = originalImage.height;

        // Iterate over 8x8 blocks in the image
        for (int y = 0; y < height; y += 8) {
            for (int x = 0; x < width; x += 8) {
                // Apply DCT to the current 8x8 block
                double[][] block = getBlock(originalImage, x, y);
                double[][] dctCoefficients = applyDCT(block);

                // Visualize the DCT coefficients in dctImage
                visualizeDCTCoefficients(dctCoefficients, dctImage, x, y);

                // Apply inverse DCT to reconstruct the block
                double[][] reconstructedBlock = applyInverseDCT(dctCoefficients);
                setBlock(reconstructedBlock, reconstructedImage, x, y);
            }
        }
    }
  
    private double[][] getBlock(RasterImage image, int startX, int startY) {
        double[][] block = new double[8][8];
       

        for (int y = startY; y < startY + 8; y++) {
            for (int x = startX; x < startX + 8; x++) {
                int pixel = image.argb[y * image.width + x];
                int grayValue = (pixel >> 16) & 0xFF;
                block[y - startY][x - startX] = grayValue - 128; // Subtract 128 for DC coefficient
            }
        }

        return block;
    }
  
    private void setBlock(double[][] block, RasterImage image, int startX, int startY) {
        for (int y = startY; y < startY + 8; y++) {
            for (int x = startX; x < startX + 8; x++) {
                double value = block[y - startY][x - startX];
                int grayValue = (int) (value + 128); // Add 128 for DC coefficient
                image.argb[y * image.width + x] = (0xFF << 24) | (grayValue << 16) | (grayValue << 8) | grayValue;
            }
        }
    }

    private double[][] applyDCT(double[][] block) {
        double[][] dctCoefficients = new double[8][8];
        double cu, cv, sum;

        for (int u = 0; u < 8; u++) {
            for (int v = 0; v < 8; v++) {
                sum = 0.0;
 
                for (int x = 0; x < 8; x++) {
                    for (int y = 0; y < 8; y++) {
                        cu = (u == 0) ? (1 / Math.sqrt(2)) : 1;
                        cv = (v == 0) ? (1 / Math.sqrt(2)) : 1;
                        sum += cu * cv * block[x][y] *
                                Math.cos(((2 * x + 1) * u * Math.PI) / 16) *
                                Math.cos(((2 * y + 1) * v * Math.PI) / 16);
                    }
                }

                dctCoefficients[u][v] = sum / 4; // Divide by 4 for normalization
            }
        }

        return dctCoefficients;
    }
  
    private double[][] applyInverseDCT(double[][] dctCoefficients) {
        double[][] block = new double[8][8];
        double cu, cv, sum;

        for (int x = 0; x < 8; x++) {
            for (int y = 0; y < 8; y++) {
                sum = 0.0;

                for (int u = 0; u < 8; u++) {
                    for (int v = 0; v < 8; v++) {
                        cu = (u == 0) ? (1 / Math.sqrt(2)) : 1;
                        cv = (v == 0) ? (1 / Math.sqrt(2)) : 1;
                        sum += cu * cv * dctCoefficients[u][v] *
                                Math.cos(((2 * x + 1) * u * Math.PI) / 16) *
                                Math.cos(((2 * y + 1) * v * Math.PI) / 16);
                    }
                }

                block[x][y] = sum / 4; // Divide by 4 for normalization
            }
        }

        return block;
    }
  
    private void visualizeDCTCoefficients(double[][] dctCoefficients, RasterImage dctImage, int startX, int startY) {
        for (int u = 0; u < 8; u++) {
            for (int v = 0; v < 8; v++) {
                double coefficient = dctCoefficients[u][v];
                int grayValue = (int) (coefficient / 8 + 128); // Scale by 8 and add 128

                int x = startX + u;
                int y = startY + v;
                dctImage.argb[y * dctImage.width + x] = (0xFF << 24) | (grayValue << 16) | (grayValue << 8) | grayValue;
            }
        }
    }
  
    public double getMSE(RasterImage originalImage, RasterImage reconstructedImage) {
        if (originalImage.width != reconstructedImage.width || originalImage.height != reconstructedImage.height) {
            throw new IllegalArgumentException("Images must have the same dimensions");
        }

        int width = originalImage.width;
        int height = originalImage.height;
        double sum = 0.0;

        for (int i = 0; i < width * height; i++) {
            int originalPixel = originalImage.argb[i];
            int reconstructedPixel = reconstructedImage.argb[i];

            int originalGray = (originalPixel >> 16) & 0xFF;
            int reconstructedGray = (reconstructedPixel >> 16) & 0xFF;

            int error = originalGray - reconstructedGray;
            sum += error * error;
        }

        return sum / (width * height);
    }
}