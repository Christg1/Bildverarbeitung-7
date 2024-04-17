// BV Ue7 SS2023 Vorgabe
//
// Copyright (C) 2023 by Klaus Jung
// All rights reserved.
// Date: 2023-03-23 		   	  	  		

package bv_ss23;
 		   	  	  		
import java.io.File;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;

import javafx.fxml.FXML;
import javafx.geometry.Point2D;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Slider;
import javafx.scene.image.ImageView;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.Region;
import javafx.stage.FileChooser;

public class DCTAppController {
	 		   	  	  		
	private static final String initialFileName = "test2.jpg";
	private static File fileOpenPath = new File(".");

	private static final DiscreteCosineTransform transform = new DiscreteCosineTransform();
	
	@FXML
	private Slider coefficientsSlider;

	@FXML
	private Label coefficientsLabel;

    @FXML
    private Slider zoomSlider;

    @FXML
    private Label zoomLabel;
 		   	  	  		
	@FXML
	private ImageView originalImageView;

	@FXML
	private ScrollPane originalScrollPane;

	@FXML
	private ImageView dctImageView;

	@FXML
	private ScrollPane dctScrollPane;

	@FXML
	private ImageView reconstructedImageView;

	@FXML
	private ScrollPane reconstructedScrollPane;

    @FXML
    private Label originalEntropyLabel;

    @FXML
    private Label dctEntropyLabel;

    @FXML
    private Label reconstructedEntropyLabel;

    @FXML
    private Label resonstructedMSELabel;

    @FXML
    private Label messageLabel;

    @FXML
    void openImage() {
    	FileChooser fileChooser = new FileChooser();
    	fileChooser.setInitialDirectory(fileOpenPath); 
    	fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("Images (*.jpg, *.png, *.gif)", "*.jpeg", "*.jpg", "*.png", "*.gif"));
    	File selectedFile = fileChooser.showOpenDialog(null);
    	if(selectedFile != null) {
    		fileOpenPath = selectedFile.getParentFile();
			RasterImage img = new RasterImage(selectedFile);
			img.convertToGray();
			img.setToView(originalImageView);
			processImages();
    		messageLabel.getScene().getWindow().sizeToScene();
    	} 		   	  	  		
    }


    @FXML
    void coefficientsChanged() {
    	processImages();
    }

    
    Point2D mousePoint;
    
    @FXML
    void mousePressed(MouseEvent event) {
    	mousePoint = new Point2D(event.getX(), event.getY());
    }
    
    @FXML
    void mouseClicked(MouseEvent event) {
    	if(Math.abs(mousePoint.getX() - event.getX()) > 5 || Math.abs(mousePoint.getY() - event.getY()) > 5) return;
    	testSelection = event.isShiftDown() ? "next" : (isTesting ? "" : "init");
    	isTesting = !isTesting || event.isShiftDown() || event.isMetaDown() || event.isAltDown() || event.isControlDown();
    	testMode = event.isMetaDown() ? "solution" : (event.isControlDown() ? "computed" : "diff");
    	processImages();
    }
    
 	@FXML
    void zoomChanged() {
    	double zoomFactor = zoomSlider.getValue();
		zoomLabel.setText(String.format("%.1f", zoomFactor));
    	zoom(originalImageView, originalScrollPane, zoomFactor);
    	zoom(dctImageView, dctScrollPane, zoomFactor);
    	zoom(reconstructedImageView, reconstructedScrollPane, zoomFactor);
    }
    
	@FXML
	public void initialize() {
		// initialize parameters
		coefficientsChanged();
		
		// load and process default image
		RasterImage img = new RasterImage(new File(initialFileName));
		img.convertToGray();
		img.setToView(originalImageView);
		processImages();
	}
	
	private void processImages() {
    	int numCoefficients = (int)coefficientsSlider.getValue();
    	coefficientsLabel.setText("" + numCoefficients);
		if(originalImageView.getImage() == null)
			return; // no image: nothing to do
		
		long startTime = System.currentTimeMillis();
		
		RasterImage origImg = new RasterImage(originalImageView);
		RasterImage dctImg = new RasterImage(origImg.width, origImg.height); 
		RasterImage recImg = new RasterImage(origImg.width, origImg.height); 
		
		transform.processDCT(origImg, dctImg, recImg, numCoefficients);

		dctImg.setToView(dctImageView);
		recImg.setToView(reconstructedImageView);
		
		String entropyFormat = "Entropy = %.3f";
		originalEntropyLabel.setText(String.format(entropyFormat, origImg.getEntropy()));
		dctEntropyLabel.setText(String.format(entropyFormat, dctImg.getEntropy()));
		reconstructedEntropyLabel.setText(String.format(entropyFormat, recImg.getEntropy()));
		
		resonstructedMSELabel.setText(String.format("MSE = %.1f", transform.getMSE(origImg, recImg)));
		
	   	messageLabel.setText("Processing time: " + (System.currentTimeMillis() - startTime) + " ms");
	   	if(isTesting)
	   		isTesting = test();
	   	else
	   		messageLabel.setEffect(null);
	}
	
	private Method testMethod = null;
	private Object testObj = null;
	private boolean isTesting = false;
	private String testSelection = "";
	private String testMode = "";

 	private boolean test() {
        try {
        	if(testMethod == null) {
        		Class<?> testClass;
        		String className = "testing.bv7c.Test";
        		try {
        			String path = System.getProperty("user.home") + File.separator + "src" + File.separator + "Java" + File.separator + "KJ_Testing.jar";
        			URL url = new File(path).toURI().toURL();
        			URLClassLoader classLoader = (URLClassLoader)ClassLoader.getSystemClassLoader();
        			Method addMethod = URLClassLoader.class.getDeclaredMethod("addURL", URL.class);
        			addMethod.setAccessible(true);
        			addMethod.invoke(classLoader, url);
        			testClass = classLoader.loadClass(className);
        		} catch (Exception e) {
            		testClass = ClassLoader.getSystemClassLoader().loadClass(className);
         		}
        		Constructor<?> constructor = testClass.getConstructor();
        		testObj = constructor.newInstance();
        		testMethod = testClass.getMethod("test", Object.class, String.class, String.class, Integer.class, String.class, String.class, String.class, String.class, String.class);
        	}
        	String image1Name = "originalImageView";
        	String image2Name = "dctImageView";
        	String image3View = "reconstructedImageView";
        	String slider1Name = "coefficientsSlider";
        	String ComboBox1Name = "";
    		testMethod.invoke(testObj, this, testSelection, testMode, 0, image1Name, image2Name, image3View, slider1Name, ComboBox1Name);
    		testSelection = "";
    		return true;
		} catch (Exception e) {
			if(testMethod != null) e.printStackTrace();
	        messageLabel.setText("No test available");
	        return false;
	    }

 	}

	private void zoom(ImageView imageView, ScrollPane scrollPane, double zoomFactor) {
		if(zoomFactor == 1) {
			scrollPane.setPrefWidth(Region.USE_COMPUTED_SIZE);
			scrollPane.setPrefHeight(Region.USE_COMPUTED_SIZE);
			imageView.setFitWidth(0);
			imageView.setFitHeight(0);
		} else {
			double paneWidth = scrollPane.getWidth();
			double paneHeight = scrollPane.getHeight();
			double imgWidth = imageView.getImage().getWidth();
			double imgHeight = imageView.getImage().getHeight();
			double lastZoomFactor = imageView.getFitWidth() <= 0 ? 1 : imageView.getFitWidth() / imgWidth;
			if(scrollPane.getPrefWidth() == Region.USE_COMPUTED_SIZE)
				scrollPane.setPrefWidth(paneWidth);
			if(scrollPane.getPrefHeight() == Region.USE_COMPUTED_SIZE)
				scrollPane.setPrefHeight(paneHeight);
			double scrollX = scrollPane.getHvalue();
			double scrollY = scrollPane.getVvalue();
			double scrollXPix = ((imgWidth * lastZoomFactor - paneWidth) * scrollX + paneWidth/2) / lastZoomFactor;
			double scrollYPix = ((imgHeight * lastZoomFactor - paneHeight) * scrollY + paneHeight/2) / lastZoomFactor;
			imageView.setFitWidth(imgWidth * zoomFactor);
			imageView.setFitHeight(imgHeight * zoomFactor);
			if(imgWidth * zoomFactor > paneWidth)
				scrollX = (scrollXPix * zoomFactor - paneWidth/2) / (imgWidth * zoomFactor - paneWidth);
			if(imgHeight * zoomFactor > paneHeight)
				scrollY = (scrollYPix * zoomFactor - paneHeight/2) / (imgHeight * zoomFactor - paneHeight);
			if(scrollX < 0) scrollX = 0;
			if(scrollX > 1) scrollX = 1;
			if(scrollY < 0) scrollY = 0;
			if(scrollY > 1) scrollY = 1;
			scrollPane.setHvalue(scrollX);
			scrollPane.setVvalue(scrollY);
		}
	}
	 		   	  	  		


 		   	  	  		
}

 		   	  	  		





