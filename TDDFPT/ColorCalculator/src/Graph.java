
import java.awt.Color;
import java.awt.Dimension;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;


import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;


/**
 * A simple demonstration application showing how to create a line chart using data from a
 * {@link CategoryDataset}.
 */
@SuppressWarnings("serial")
public class Graph extends ApplicationFrame{

    /**
     * Creates a new demo.
     *
     * @param title  the frame title.
     * @throws IOException 
     */
    public Graph(  String title) throws IOException {
        super(title);
    	
          XYDataset dataset = createDataset("color.dat");
          JFreeChart chart = createChart(dataset);
          ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(500, 270));
        setContentPane(chartPanel);
    }
    public Graph(String title,String fileName) throws IOException {
       super(title);
    	
          XYDataset dataset = createDataset(fileName);
          JFreeChart chart = createChart(dataset);
          ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(500, 270));
       setContentPane(chartPanel);
         
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
       setVisible(true);

    }
    public Graph(  String title, Double[][]arrayData) throws IOException {
        super(title);
          XYDataset dataset = createDataset(arrayData);
          JFreeChart chart = createChart(dataset);
          ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        setContentPane(chartPanel);
         
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);

    }
    
    

    /**
     * Creates a sample dataset.
     * 
     * @return The dataset.
     */
    
    public XYDataset createDataset(String dataFileName) throws IOException {
    	
    	XYSeries series1 = new XYSeries("First"); 
    	
    	XYSeries series2 = new XYSeries("Second");
    	
    	//DefaultCategoryDataset dataset = new DefaultCategoryDataset();
    	byte[] buffer = new byte[(int) new File(dataFileName).length()];
        BufferedInputStream f = null;
        try {
                f = new BufferedInputStream(new FileInputStream(dataFileName));
                f.read(buffer);
        } catch (FileNotFoundException e) {
                e.printStackTrace();
        } catch (IOException e) {
                e.printStackTrace();
        }  finally {
                if (f != null) try { f.close(); } catch (IOException ignored) { }
        }

        String data = new String(buffer);
        
        data = data.replace("\t"," ");

        String[] lambda_val = data.split("\n");

        Double[][] colorData = new Double[lambda_val.length][2];


        try {
                for (int i=0; i<lambda_val.length; i++){
                        String[] lambda_color = lambda_val[i].split(" ");
                        colorData[i][0] = Double.parseDouble(lambda_color[0]); 
                        colorData[i][1] = Double.parseDouble(lambda_color[1]);
                        series2.add(colorData[i][0], colorData[i][1]);
                }

        } catch (Exception e){
       	 System.err.println("FILENAME: "+dataFileName);
                e.printStackTrace();
        }
//    	       BufferedReader bReader =new BufferedReader(new FileReader(dataFileName));
//    	          String s;
//          while (!(s=bReader.readLine()).equals("") && !s.equals(null)){
//        	          	  
//        	  s = s.replace("\t"," ");                      
//        	  String datavalue [] = s.split(" ");
//    	              String category = datavalue[0];
//    	              System.out.println(category);
//    	              String value = datavalue [1];
//    		              dataset.addValue(Double.parseDouble(value), series1, category);
//    		if (datavalue[0]=="780") break;            
//          }
//    	            bReader.close();
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series1);
        dataset.addSeries(series2);
    		return dataset;
    	    }
    
    public XYDataset createDataset(Double[][] dataArray)  {
    	XYSeries series1 = new XYSeries("First");    
    	
    	XYSeries series2 = new XYSeries("Second");
    	
    	
    	//DefaultCategoryDataset dataset=new DefaultCategoryDataset();
		for(int i =0; i<dataArray.length;i++){
			series2.add(dataArray[i][0], dataArray[i][1]);
    		            }
		XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series1);
        dataset.addSeries(series2);            
    		        return dataset;
    	    }
    
    /**
     * Creates a sample chart.
     * 
     * @param dataset  a dataset.
     * 
     * @return The chart.
     */
 private JFreeChart createChart(  XYDataset dataset) {
        
        // create the chart...
          JFreeChart chart = ChartFactory.createXYLineChart(
            "Line Chart Demo 6",      // chart title
            "X",                      // x axis label
            "Y",                      // y axis label
            dataset,                  // data
            PlotOrientation.VERTICAL,
            true,                     // include legend
            true,                     // tooltips
            false                     // urls
        );

        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
        chart.setBackgroundPaint(Color.white);

//          StandardLegend legend = (StandardLegend) chart.getLegend();
  //      legend.setDisplaySeriesShapes(true);
        
        // get a reference to the plot for further customisation...
          XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.lightGray);
    //    plot.setAxisOffset(new Spacer(Spacer.ABSOLUTE, 5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        
          XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, false);
        renderer.setSeriesShapesVisible(1, false);
        plot.setRenderer(renderer);

        // change the auto tick unit selection to integer units only...
          NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        // OPTIONAL CUSTOMISATION COMPLETED.
                
        return chart;
        
    }
    
    /**
     * Starting point for the demonstration application.
     *
     * @param args  ignored.
     * @throws IOException 
     */
    public static void main(  String[] args) throws IOException {

    	Graph demo = new Graph("Line Chart Demo 6");
        demo.pack();
        RefineryUtilities.centerFrameOnScreen(demo);
        demo.setVisible(true);

    }

}
