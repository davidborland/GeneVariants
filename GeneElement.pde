/*=========================================================================
 
 Name:        GeneElement.pde
 
 Author:      David Borland, The Renaissance Computing Institute (RENCI)
 
 Copyright:   The Renaissance Computing Institute (RENCI)
 
 Description: Base classes for elements located at some position along a gene.
 
 =========================================================================*/
 
 
 abstract class GeneElement {  
  // Start, end, and length
  int start;
  int end;
  int length;
  
  
  // Constructors
  GeneElement(int _start, int _end) {
    start = _start;
    end = _end;
    
    ComputeLength();
    
    // XXX: Gene element should take the Geometry as an argument, and determine what its transcript->gene offset is...
  }
  
  GeneElement() {
    start = 1;
    end = 1;
    
    ComputeLength();
  }
  
  
  // Draw the element
  abstract void Draw();
  
  
  // Compute the length
  protected void ComputeLength() {
    // Length is inclusive, as we are dealing with discrete positions
    length = end - start + 1;
  }
}


abstract class GeneElements {   
  // Array list of elements
  ArrayList elements;
  
  
  // Constructor
  GeneElements() {
    // Create array list
    elements = new ArrayList();  
  }
  
  
  // Draw the elements
  void Draw() {       
    for (int i = 0; i < elements.size(); i++) {
      GeneElement e = (GeneElement) elements.get(i);
      e.Draw();
    }
  }
  
  
  // Add an element 
  protected void AddElement(String values[], HashMap<String, Integer> columnMap) {
  }
  
  
  // Read the data from a results text file
  void ReadFile(String fileName) {
      // Load the file into an array of strings
    String fileLines[] = loadStrings(fileName);
    
    
    // Create a hashmap so we can index by column name
    HashMap<String, Integer> columnMap = new HashMap();
    String columnNames[] = split(fileLines[0], '\t');
    for (int i = 0; i < columnNames.length; i++) {
      columnMap.put(columnNames[i], new Integer(i));
    }
        
    // Loop over the lines in the file, skipping the first
    for (int i = 1; i < fileLines.length; i++) {
      // Split the current line by tabs
      String values[] = split(fileLines[i], '\t');      
      
      // Check for a comment
      if (values[0].charAt(0) == '#') continue;
      
      AddElement(values, columnMap);
    }
  }
}
