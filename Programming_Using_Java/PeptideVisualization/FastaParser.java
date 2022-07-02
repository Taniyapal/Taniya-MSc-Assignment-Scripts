package PeptideVisualization;

import java.io.File;
import java.io.FileNotFoundException;

import java.io.IOException;
import java.util.HashMap;
import java.util.Scanner;



public class FastaParser {
    
    String fileName;
    String header;
    String sequence;
    
    String line;
    String accesion;
    HashMap <String, String> h = new HashMap <String,String>();
    
    
    public FastaParser(String file) {
        fileName = file;
        File fastafile= new File (fileName);
        if(!fastafile.exists()) {
            System.out.println("Fasta file not found!");
        }    
        
    }
    
    public void FastaLineReader(){
        Scanner myscanner;
        File fastafile= new File(fileName);
        try {
            myscanner= new Scanner(fastafile);
            while(myscanner.hasNextLine()){
                String line= myscanner.nextLine();
                if (line.startsWith(">")){
                    header=line;
                    String [] array=header.split("\\|");
                    accesion=array[1];
                    
                    sequence="";
                } else {
                    sequence=sequence+line;
                }
                h.put(header, sequence);
            }
            myscanner.close();
        } catch (FileNotFoundException ex) {
            System.out.println(ex);
        }
        
        
        
    }
    
    
}
