import javax.swing.JList;

import java.io.File;
import java.util.ArrayList;


class GUITools {
	
	
	public static JList getFolderList(String folderPath){
		
		if(folderPath==null) folderPath=".";
		File folder = new File(folderPath);
	    File[] listOfFiles = folder.listFiles();

	    ArrayList<String> datfiles = new ArrayList<String>();
	    for (int i = 0; i < listOfFiles.length; i++)
	      if (listOfFiles[i].isFile()) 
	        if(listOfFiles[i].getName().endsWith(".dat"))
	        	datfiles.add(listOfFiles[i].getName());
		
		return new JList(datfiles.toArray(new String[datfiles.size()]));
	}

}
