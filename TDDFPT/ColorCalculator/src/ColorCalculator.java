import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.swing.JPanel;


public class ColorCalculator {

	private static final double[][] XYZtoRGBMatrix = {
		{3.240479, -1.537150, -0.498535},
		{-0.969256, 1.875992, 0.041556},
		{0.055648, -0.204043, 1.057311}
	};

	//	private static double AMOUNT = 0.95;

	private double[][] light=null;
	private double[][] color=null;
	double[][] CIEx = null;
	double[][] CIEy = null;
	double[][] CIEz = null;
	double[][] result = null;

	public ColorCalculator(){
		CIEx = parseColor("CIEx");
		CIEy = parseColor("CIEy");
		CIEz = parseColor("CIEz");
		loadLight();
	}

	public void loadLight(){
		light=parseColor("D65.light");
	}

	public void loadAbsorbSpectrum(String filename){
		color = parseColor(filename);    
        //        System.out.println("begin interpolate");
                color = interpolateColor(color);	// By XC.G: in order to adapt color files with different grids 
	}

//  By XC.G  18.may.2012
        public static double[][] interpolateColor(double[][] color){
        	// first step :  sort, ascending
		for (int i = 0; i < color.length-1;i++){
			int min=i;
			for(int j = i+1; j < color.length; j++){
				if(color[j][0] < color[min][0]) min = j;
			}
			double temp0 = color[min][0], temp1 = color[min][1];
                        color[min][0]=color[i][0];
			color[min][1]=color[i][1];
			color[i][0]=temp0;
			color[i][1]=temp1;
		}

                // System.out.println("end sorting");
		// sencond step : interpolate
		double[][] color_inter = new double[81][2];	// from 380 to 780, step 5, so the number of grid is 81
                int if_end=0;			// flag to judge if we already got all the points for color_inter
		int i_wl=0;			// color_inter[i_wl][0]=380+5*i_wl
		for (int i = 0; i < 81; i++) color_inter[i][0]=380+i*5;
		double[] head = {0,0};
		double[] tail = {0,0};		// head and tail is a train moving on color, head = color[i] and tail = color[i-1]; when color_inter[j][0] 
						// is between the head[0] and tail[0], do the interpolation
		for(int i = 0; i < color.length-1;i++){
			tail = head;
                        head = color[i];
			while(head[0]>=color_inter[i_wl][0]&&tail[0]<color_inter[i_wl][0]){	// color_inter[i_wl] is in the train
				color_inter[i_wl][1]=((head[0]-color_inter[i_wl][0])*tail[1]+(color_inter[i_wl][0]-tail[0])*head[1])/(head[0]-tail[0]);
                                System.out.println(i_wl);
                                System.out.println(i);
                                System.out.println(color_inter[i_wl][0]);
                                System.out.println(head[0]);
                                System.out.println(tail[0]);
				i_wl=i_wl+1;
				if(i_wl>=81) { if_end=1; break;}
			}
                        if(if_end==1) break;
		}
		return color_inter;
        }
// End of XC.G

	public static double[][] parseColor(String filename){

		//read file from input
		byte[] buffer = new byte[(int) new File(filename).length()];
		BufferedInputStream f = null;
		try {
			f = new BufferedInputStream(new FileInputStream(filename));
			f.read(buffer);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (f != null) try { f.close(); } catch (IOException ignored) { }
		}

		String data = new String(buffer);

		data = data.replace("\t"," ");

		String[] lambda_val = data.split("\n");

		double[][] colorData = new double[lambda_val.length][2];


		try {
			for (int i=0; i<lambda_val.length; i++){
				String[] lambda_color = lambda_val[i].split(" ");
				colorData[i][0] = Double.parseDouble(lambda_color[0]); 
				colorData[i][1] = Double.parseDouble(lambda_color[1]); 
			}

		} catch (Exception e){
			System.err.println("ERROR READING "+filename);
			e.printStackTrace();
		}


		return colorData;
	}


	private static double findMax(double[][] data){
		double max = 0;

		for (int i=0; i<data.length; i++)
			if (data[i][1]>max)
				max = data[i][1];
		return max;
	}


	static double[] XYZtoRGB(double X, double Y, double Z, double[][] m){
		double[] rgb = new double[3];

		rgb[0] = m[0][0]*X + m[0][1]*Y + m[0][2]*Z;
		rgb[1] = m[1][0]*X + m[1][1]*Y + m[1][2]*Z;
		rgb[2] = m[2][0]*X + m[2][1]*Y + m[2][2]*Z;

		return rgb;
	}

	//	public static void main(String[] args){
	//
	//		double[][] color_dat = parseColor("color.dat");
	//		double[][] Io = parseColor("WhiteLight.dat");
	//		double[][] CIEx = parseColor("CIEx.dat");
	//		double[][] CIEy = parseColor("CIEy.dat");
	//		double[][] CIEz = parseColor("CIEz.dat");
	//		double d=AMOUNT/findMax(color_dat);
	//		double[][] I=calcIntensitySpec(Io,color_dat,d);
	//		double X=calcXYZ(I,CIEx);
	//		double Y=calcXYZ(I,CIEy);
	//		double Z=calcXYZ(I,CIEz);
	//		double Ylum = calcXYZ(Io, CIEy);
	//		//         Ylum *= (1-AMOUNT);
	//
	//		System.out.println("X= "+ X+ ", Y= " +Y+", Z= " +Z+", Ylum= "+Ylum);
	//		X=X/Ylum;
	//		Y=Y/Ylum;
	//		Z=Z/Ylum;
	//		//         Ylum/=Ylum;
	//		double tot = X+Y+Z;
	//		double Xnorm = X / tot, Ynorm = Y / tot, Znorm = Z/ tot;
	//		//         for (int i=0; i<color_dat.length; i++)
	//		//                 System.out.println(color_dat[i][0]+"\t"+color_dat[i][1]);
	//		//         System.out.println("--------------------------------------------");
	//		//         for (int i=0; i<Io.length; i++) 		
	//		//         		 System.out.println(Io[i][0]+"\t"+Io[i][1]);
	//		//         System.out.println("--------------------------------------------");
	//		//         for (int i=0; i<x.length; i++) 		
	//		//         		 System.out.println(x[i][0]+"\t"+x[i][1]);
	//		//         for (int i=0; i<I.length; i++) 		
	//		//     		 System.out.println(I[i][0]+"\t"+I[i][1]);
	//		System.out.println("X= "+ X+ ", Y= " +Y+", Z= " +Z);
	//		System.out.println("Xn= "+ Xnorm+ ", Yn= " +Ynorm+", Zn= " +Znorm+", Ylum= "+Ylum);
	//		//         double[] rgb = XYZtoRGB(Xnorm, Ynorm, Znorm, XYZtoRGBMatrix);
	//		double[] rgb = XYZtoRGB(Xnorm, Ynorm, Znorm, XYZtoRGBMatrix);
	//
	//		System.out.println("R= "+rgb[0]+", G= "+rgb[1]+", B= "+rgb[2]);
	//
	//		//APPLY GAMMA
	//		for(int j=0; j<3; j++){
	//			if (rgb[j]>1) rgb[j]=1;
	//
	//			double a = 0.055;
	//			if (rgb[j]<0.0031308)
	//				rgb[j]*=12.92;
	//			else
	//				rgb[j] = (1+a) * Math.pow(rgb[j], 1/2.4) - a; 
	//		}
	//
	//		System.out.println("R= "+rgb[0]+", G= "+rgb[1]+", B= "+rgb[2]);
	//		Color color = new Color((float)(rgb[0]), (float)(rgb[1]), (float)(rgb[2]));
	//		JFrame frame = new JFrame();
	//		JPanel panel = new JPanel();
	//		panel.setBackground(color);
	//		frame.setMinimumSize(new Dimension(600,300));
	//		frame.getContentPane().add(panel);
	//		frame.pack();
	//		frame.setVisible(true);
	//	}

	public static double[][] calcIntensitySpec(double[][] Io, double[][] c, double x){
		double[][] I = new double[Io.length][2];
		double cmax = findMax(c);
		for(int i=0;i<Io.length;i++){
			I[i][0]=Io[i][0];
			//		 I[i][1]=Io[i][1]*(1-c[i][1]*dist);
			
//			I[i][1]=Io[i][1]*(1-c[i][1]*x);
			I[i][1]=Io[i][1]*Math.exp( c[i][1]*x/cmax );
		}	 

		return I;
	}


	public static double calcXYZ(double[][] I, double[][] CIE){
		double W = 0;

		for(int i=0;i<I.length;i++){
			W+= I[i][1]*CIE[i][1];
		}	 

		return 5*W;
	}


	//	public static double[] calculateColor(double amount, double brightnessOffset){
	//
	//		double[][] color_dat = parseColor("color.dat");
	//		double[][] Io = parseColor("WhiteLight.dat");
	//		double[][] CIEx = parseColor("CIEx.dat");
	//		double[][] CIEy = parseColor("CIEy.dat");
	//		double[][] CIEz = parseColor("CIEz.dat");
	//		double d=amount/findMax(color_dat);
	//		double[][] I=calcIntensitySpec(Io,color_dat,d);
	//		double X=calcXYZ(I,CIEx);
	//		double Y=calcXYZ(I,CIEy);
	//		double Z=calcXYZ(I,CIEz);
	//		double Ylum = calcXYZ(Io, CIEy);
	//		Ylum *= (1-AMOUNT);
	//		System.out.println("X= "+ X+ ", Y= " +Y+", Z= " +Z+", Ylum= "+Ylum);
	//		//     X=X/Ylum;
	//		//     Y=Y/Ylum;
	//		//     Z=Z/Ylum;
	//		//     Ylum/=Ylum;
	//		double tot = X+Y+Z;
	//		double Xnorm = X / tot, Ynorm = Y / tot, Znorm = Z/ tot;
	//		//     for (int i=0; i<color_dat.length; i++)
	//		//             System.out.println(color_dat[i][0]+"\t"+color_dat[i][1]);
	//		//     System.out.println("--------------------------------------------");
	//		//     for (int i=0; i<Io.length; i++) 		
	//		//     		 System.out.println(Io[i][0]+"\t"+Io[i][1]);
	//		//     System.out.println("--------------------------------------------");
	//		//     for (int i=0; i<x.length; i++) 		
	//		//     		 System.out.println(x[i][0]+"\t"+x[i][1]);
	//		//     for (int i=0; i<I.length; i++) 		
	//		// 		 System.out.println(I[i][0]+"\t"+I[i][1]);
	//		System.out.println("X= "+ X+ ", Y= " +Y+", Z= " +Z);
	//		Ynorm*= (1+brightnessOffset);
	//		System.out.println("Xn= "+ Xnorm+ ", Yn= " +Ynorm+", Zn= " +Znorm+", Ylum= "+Ylum);
	//		//     double[] rgb = XYZtoRGB(Xnorm, Ynorm, Znorm, XYZtoRGBMatrix);
	//		double[] rgb = XYZtoRGB(Xnorm, Ynorm, Znorm, XYZtoRGBMatrix);
	//
	//		System.out.println("R= "+rgb[0]+", G= "+rgb[1]+", B= "+rgb[2]);
	//
	//		return rgb;
	//	}


	public double[] calculateIlluminatedColor(double amount, double brightnessOffset){

		if (color==null) {
			double[] white = new double[3];
			white[0]=1;
			white[1]=1;
			white[2]=1;
			return white;
		}

		double[][] Io = light;

//		double d=amount/findMax(color);
		if (amount==0) amount=0.00001;
		double d=Math.log(amount);
		double[][] I=calcIntensitySpec(Io,color,d);
		result = I;
		double X=calcXYZ(I,CIEx);
		double Y=calcXYZ(I,CIEy);
		double Z=calcXYZ(I,CIEz);
		double Ylum = calcXYZ(Io, CIEy);
		//     Ylum *= (1-AMOUNT);
		//		System.out.println("X= "+ X+ ", Y= " +Y+", Z= " +Z+", Ylum= "+Ylum);
		X=X/Ylum;
		Y=Y/Ylum;
		Z=Z/Ylum;
		//     Ylum/=Ylum;
		//double tot = X+Y+Z;

		double[] rgb = XYZtoRGB(X, Y, Z, XYZtoRGBMatrix);

		System.out.println("R= "+rgb[0]+", G= "+rgb[1]+", B= "+rgb[2]);

		return rgb;
	}


	public static double[][] createGaussian(double lambda1, double width1, double lambda2, double width2){
		double[][] Io = parseColor("D65.light");
		//		Double lambda1=680.0;//red
		//		Double lambda2=450.0;//blue
		//Double lambda2=530.0;//green

		double[][] out = new double[Io.length][2];
		for(int i=0;i<Io.length;i++){
			out[i][0] = Io[i][0];
			out[i][1]=(Math.exp(-Math.pow(((Io[i][0]-lambda1)/width1),2.0))+
					Math.exp(-Math.pow(((Io[i][0]-lambda2)/width2),2.0)));
		}

		return out;
	}

	public JPanel createAbsorbptionGraph(){
		JPanel graph = null;
		if(color!=null){
			GraphFactory gf = new GraphFactory("Absorbption Spectrum");
			gf.addDataSet("Absorbption spectrum", color);
			graph = gf.getGraph();
		} else {
			graph = new JPanel();
		}

		return graph;
	}
	
	public JPanel createResultLightGraph(){
		JPanel graph = null;
		if(result!=null){
			GraphFactory gf = new GraphFactory("Result light Spectrum");
			gf.addDataSet("visible spectrum", result);
			graph = gf.getGraph();
		} else {
			graph = new JPanel();
		}

		return graph;
	}
	
	public JPanel createGraph(){
		JPanel graph = null;
		
		if(result!=null && color!=null){
			GraphFactory gf = new GraphFactory("");
			gf.addDataSet("Absorbption spectrum", color);
			gf.addDataSet("Visible spectrum", result);
			graph = gf.getGraph();
		} else {
			graph = new JPanel();
		}

		return graph;
	}



	private static Double[][] toDouble(double[][] data){
		Double[][] out = new Double[data.length][];

		for (int i=0; i<data.length; i++){
			out[i] = new Double[data[i].length];
			for (int k=0; k<data[i].length; k++)
				out[i][k]= new Double(data[i][k]);
		}

		return out;
	}

}
