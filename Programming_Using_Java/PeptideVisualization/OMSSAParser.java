package PeptideVisualization;

import au.com.bytecode.opencsv.CSVReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


import java.util.Arrays;

import java.util.List;






public class OMSSAParser {
    
    String filename;
    
    List body;
    
    
    //Vector rows;
   
    
    
    public OMSSAParser(String file) {
        filename = file;
        File OMSSAfile= new File (filename);
        if(!OMSSAfile.exists()) {
            System.out.println("OMSSA file not found!");
        }    
        
    }
    
    public void OMSSAFileReader() throws FileNotFoundException, IOException  {
        
        try {
            File OMSSAfile= new File (filename);
            CSVReader reader = new CSVReader (new FileReader(OMSSAfile));
            body=reader.readAll();
            
            
       
        } catch (Exception e) { 
            System.out.println("file not found"+ e); 
 
        } 
        
            
        
        
        
 
            
            
        
        
            
    }

}        
            
            
            
        
        
    

