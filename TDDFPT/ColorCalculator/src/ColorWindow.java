import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextArea;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;


public class ColorWindow {

	private double[] rgb = new double[3];
	
	private JFrame frame;
	private JPanel colorPanel;
	private JTextArea rgbvalues;
	private JTextArea rgbhexvalues;
	private JSlider amountBar;
	private JSlider brightnessBar;
	private Box graphBox;
	
	private double AMOUNT=0.5, BRIGHTNESS=0;
	private Boolean GAMMA=false, MERGE=false;
	private String SPECTRUMFILE="", PREV_SPECTRUMFILE="";
	private ColorCalculator cc = null;

	public ColorWindow(){

		initComponents();
		
		cc=new ColorCalculator();
	}
	
	
	public static void main(String args[]){
		new ColorWindow();
	}


	private double applyGamma(double c){
		if (c>1) c=1;

		double a = 0.055;
		if (c<0.0031308)
			c*=12.92;
		else
			c = (1+a) * Math.pow(c, 1/2.4) - a; 
		
		if (c>1) c=1;
		
		return c;
	}
	

	private void initComponents(){
		frame = new JFrame();
		frame.setLayout(new BorderLayout());
		
		JCheckBox gammayesno = new JCheckBox("Use Gamma Correction" , GAMMA);
		gammayesno.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
					GAMMA = ((JCheckBox)e.getSource()).isSelected();
					action_calculateColor();
				}
			});
		
		JCheckBox merge = new JCheckBox("Merge graphs", MERGE);
		merge.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
					MERGE = ((JCheckBox)e.getSource()).isSelected();
					drawGraphs(true);
				}
			});
		
		Box toolsBox = Box.createHorizontalBox();
		toolsBox.add(gammayesno);
		toolsBox.add(merge);
		
		
		amountBar = new JSlider(0, 1000, 500);
		amountBar.addChangeListener(
				new ChangeListener() {
					@Override
					public void stateChanged(ChangeEvent e) {
						JSlider source = (JSlider)e.getSource();
						if (!source.getValueIsAdjusting())
							AMOUNT = (double)(source.getModel().getValue())/1000;
							action_calculateColor();
					}
				});
		amountBar.setPaintLabels(true);
		
		brightnessBar = new JSlider(JSlider.VERTICAL, -100,100,0);
		brightnessBar.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				JSlider source = (JSlider)e.getSource();
				if (!source.getValueIsAdjusting())
					BRIGHTNESS = (double)(source.getModel().getValue())/1000;
					action_calculateColor();
			}
		});
		
		
		JList folderList = GUITools.getFolderList(null);
		JPanel listPanel = new JPanel();
//		listPanel.setPreferredSize(new Dimension(120,500));
		listPanel.add(folderList);
		folderList.addListSelectionListener(
				new ListSelectionListener() {
					@Override
					public void valueChanged(ListSelectionEvent e) {
						PREV_SPECTRUMFILE = SPECTRUMFILE;
						SPECTRUMFILE = ((JList)e.getSource()).getSelectedValue().toString();
						action_calculateColor();
					}
				});
		JScrollPane scroll = new JScrollPane(listPanel);
		scroll.setAutoscrolls(true);
		scroll.setPreferredSize(new Dimension(125,500));
		
		Box amBox = Box.createHorizontalBox();
		amBox.add(new JLabel("Amount"));
		amBox.add(amountBar);
		
		Box rightBox = Box.createVerticalBox();
		rgbvalues = new JTextArea();
		rgbhexvalues = new JTextArea();
		
		rightBox.add(rgbvalues);
		rightBox.add(rgbhexvalues);
		colorPanel = new JPanel();
		colorPanel.setMinimumSize(new Dimension(100,100));
		colorPanel.setPreferredSize(new Dimension(150,400));
		rightBox.add(colorPanel);
		
		graphBox = Box.createVerticalBox();
		
		frame.setMinimumSize(new Dimension(600,400));
		frame.getContentPane().add(rightBox, BorderLayout.EAST);
		frame.getContentPane().add(toolsBox, BorderLayout.NORTH);
		frame.getContentPane().add(graphBox, BorderLayout.CENTER);
		frame.getContentPane().add(amBox, BorderLayout.SOUTH);
//		frame.getContentPane().add(brightnessBar, BorderLayout.WEST);
//		frame.getContentPane().add(listPanel, BorderLayout.WEST);
		frame.getContentPane().add(scroll, BorderLayout.WEST);
		
		frame.setTitle("Color viewer version 5");
		
		frame.pack();
		frame.setVisible(true);
		
		
		
	}
	
	
	private void action_calculateColor(){
		
		if(!SPECTRUMFILE.equals(PREV_SPECTRUMFILE)){
			cc.loadAbsorbSpectrum(SPECTRUMFILE);
		}
		
		rgb = cc.calculateIlluminatedColor(AMOUNT, BRIGHTNESS);
		
		System.out.println("Amount: "+AMOUNT+" , BRIGHTNESS: "+BRIGHTNESS+" GAMMA correction: "+GAMMA);
		
		if(GAMMA)
			for (int i=0; i<3; i++)
				rgb[i] = applyGamma(rgb[i]);
		
		// when transforming from XYZ to RGB some values can be negative or greater than 1
		// this means the color cannot be exactly represented in RGB space
		for (int i=0;i<rgb.length;i++){
			 if(rgb[i]<0)
				 rgb[i]=0.0;
			 if(rgb[i]>1)
				 rgb[i]=1.0;
		 }
		
		Color color = new Color ((float)rgb[0], (float)rgb[1], (float)rgb[2]);
		colorPanel.setBackground(color);
		rgbvalues.setText("RGB: "+ color.getRed() + " " + color.getGreen() + " " + color.getBlue());
		rgbhexvalues.setText("RGB: #"+ Integer.toHexString(color.getRed()).toUpperCase() + Integer.toHexString(color.getGreen()).toUpperCase() + Integer.toHexString(color.getBlue()).toUpperCase() );
		
		drawGraphs(false);
		
		frame.validate();
	}
	
	
	private void drawGraphs(Boolean validate){
		graphBox.removeAll();
		if (MERGE){
			graphBox.add(cc.createGraph());
		} else {
			graphBox.add(cc.createAbsorbptionGraph());
			graphBox.add(cc.createResultLightGraph());
		}
		
		if(validate)
			frame.validate();
	}
	

}
