/*=========================================================================
 
 Name:        GeneVariants.pde
 
 Author:      David Borland, The Renaissance Computing Institute (RENCI)
 
 Copyright:   The Renaissance Computing Institute (RENCI)
 
 Description: Processing sketch to visualize gene-variant data.  Based on
              plotit.py from Chris Bizon.
 
 =========================================================================*/
 
 
 import java.awt.event.*;
 
 
// Layout
int border = 20;  // border around screen in pixels
float sceneScale = 1.0;

float intronScale = 0.0;


// Objects
Geography geography;

Sequences sequences;
Regions regions;
Results results;
Hotspots hotspots;


// Text
PFont font;


void setup() {
  // Processing setup
  size(1024, 512);
  
  colorMode(RGB, 1.0);
  
  font = createFont("Arial", 16, true);
  
  addMouseWheelListener(new MouseWheelListener() { 
    public void mouseWheelMoved(MouseWheelEvent mwe) { 
      mouseWheel(mwe.getWheelRotation());
  }}); 
  
  // Read geography data
  geography = new Geography();
  geography.ReadFile("geographyGRN.txt");

  // Read sequence data
  sequences = new Sequences();
//  sequences.SetGeography(geography);
  sequences.ReadFile("geographyGRN.txt");
  
  // Read regions data
  regions = new Regions();
  regions.ReadFile("regionsGRN.txt");
  
  // Read results data
  results = new Results();
  results.ReadFile("results_GRN.txt");
  
  // Read hotspots data
  hotspots = new Hotspots();
  hotspots.ReadFile("hotspots_GRN.txt");
}


void mouseWheel(int delta) {
  intronScale = constrain(intronScale + max(0.01, intronScale) * 0.1 * -delta, 0.0, 1.0);
}


void draw() {
// Just in case we get stuck in a weird off-screen location
//  frame.setLocation(0, 0);

  // General setup  
  background(0.95, 0.95, 0.9);
  smooth();


  // Set projection based on length of transcript
  setProjection((int)geography.GetScreenWidth(intronScale), border);
//setProjection(sequences.codingLength + 50 , border);

  // Draw the objects
  sequences.Draw();
  regions.Draw();
  results.Draw();
  hotspots.Draw();
  
  
  // User interaction 
  int x = constrain((int)screen2transcript((float)mouseX), sequences.transcriptStart, sequences.transcriptEnd);

  stroke(0.0, 0.0, 0.0, 0.5);
  strokeWeight(1);
  
  line(x, height / 2 - 20, x, height / 2 + 180);
  
    
  resetMatrix();
  
  textFont(font);
  fill(0);  
  
  String s = "Position";
  s += "\nGenomic: " + geography.GeneToGenomic(x);
  s += "\nGene: " + x;
  s += "\nTranscript: " + geography.GeneToTranscript(x);
  s += "\nCoding: " + geography.TranscriptToCoding(geography.GeneToTranscript(x));
//  text("Position:\n\tGenomic: " + geography.GeneToTranscript(x), Gene position: " + x + "\nTranscript Position: " + geography.GeneToTranscript(x), border, border * 2);
text(s, border, border * 2);
}


void setProjection(int transcriptLength, int border) {
  // Translate the number of pixels indicated for the border
  translate(border, 0);

  // Scale based on the transcript length and border
  sceneScale = (float)(width - border * 2) / transcriptLength;
  scale(sceneScale, 1.0);
}


float screen2transcript(float xScreen) {
  return (xScreen - border) / sceneScale;
}


class Sequences extends GeneElements {
  // Start and end positions for the transcript, aka pre-mRNA
  int transcriptStart;
  int transcriptEnd;
  int transcriptLength;

  // Start and end positions for the coding region
  int codingStart;
  int codingEnd;  
  int codingLength;


  // Constructor
  Sequences() {
    super();
    
    // Just set some default values
    transcriptStart = 1;
    transcriptEnd = 1;
    
    // Setting length to the difference plus one, as these are discrete positions, 
    // so for example a transcript that begins at 4 and ends at 7 would include positions
    // 4, 5, 6, and 7, therefore having a length of 7 - 4 + 1 = 4.
    transcriptLength = transcriptEnd - transcriptStart + 1;

    codingStart = 1;
    codingEnd = 1;
    
    // See transcript length above
    codingLength = codingEnd - codingStart + 1;  
  }


  // Draw the geography
/*
  void Draw() {    
    for (int i = 0; i < sequences.size(); i++) {
      Sequence s = (Sequence) sequences.get(i);
      s.Draw();
    }
*/
/*
    fill(0.0, 0.0, 1.0);
    
    float hTranscript = 2.0;
    float hCoding = 10.0;
    
    rect(transcriptStart, height / 2.0 - hTranscript / 2.0, transcriptLength, hTranscript);
    rect(codingStart, height / 2.0 - hCoding / 2.0, codingLength, hCoding);
*/
//  }


  // Read the data from a geography text file
  void ReadFile(String fileName) {
    // Load the file into an array of strings
    String lines[] = loadStrings(fileName);


    // Looks like the coding start and end are always the two values on line 6 in the file.
    // Also looks like tabs are always the delimiter.
    String codingLine[] = split(lines[5], '\t');

    // Do a simple test until we have a better specification
    if (codingLine.length == 2) {      
      codingStart = int(codingLine[0]);
      codingEnd = int(codingLine[1]);
      codingLength = codingEnd - codingStart + 1;
    }
    else {
      println("Geography::ReadFile() : Coding start and end not where expected");
      println(codingLine);

      return;
    } 


    // Looks like the transcript start is always 1
    transcriptStart = 1;

    // Looks like the transcript end is the last number on the last line
    String transcriptLine[] = split(lines[lines.length - 1], '\t');
    transcriptEnd = int(transcriptLine[4]);

    transcriptLength = transcriptEnd - transcriptStart + 1;
    
    
String geneStartLine[] = split(lines[6], '\t'); 
int geneOffset = int(geneStartLine[0]) - 1;
    
    // Find the sequences
    for (int i = 6; i < lines.length; i++) {
      // Split the current line by tabs
      String sequenceLine[] = split(lines[i], '\t');
      
      // Get the start and end positions of this sequence.
//      int sequenceStart = int(sequenceLine[3]);
//      int sequenceEnd = int(sequenceLine[4]);

      int sequenceStart = int(sequenceLine[0]) - geneOffset;
      int sequenceEnd = int(sequenceLine[1]) - geneOffset;
      
      println(sequenceStart);
      println(sequenceEnd);
      println();
      
      // Add a sequence 
      elements.add(new Sequence(sequenceStart, sequenceEnd, sequenceLine[2]));
    }
    
    Sequence s = (Sequence) elements.get(elements.size() - 1);
    transcriptEnd = s.end;
    
    transcriptLength = transcriptEnd - transcriptStart + 1;
  }
}


class Sequence extends GeneElement {  
  // Sequence type, either UTR, intron, or exon
  // XXX: Should probably make this an enum
  String type;
  
  
  // Constructor
  Sequence(int _start, int _end, String _type) {
    super(_start, _end);
    
    type = _type;
  }
 
 
  // Draw the sequence
  void Draw() {    
    float h = 0.0;    
    
    if (type.contains("UTR")) { 
      noStroke();
      fill(0.0, 0.0, 1.0, 0.5);
      
      h = 2.0;      
    }
    else if (type.equals("exon")) {
      noStroke();
      fill(0.0, 0.0, 1.0, 1.0);
    
      h = 16.0;
    }
    else if (type.equals("intron")) {
      stroke(0.0, 0.0, 1.0, 0.75);
      strokeWeight(1);
      noFill();
      
      h = 8.0;
    }
      
//    rect(start, height / 2.0 - h / 2.0, length, h);
double x = geography.GeneToScreen(start, intronScale);

double l = length;
if (type.equals("intron")) l = length * intronScale;

rect((int)x, height / 2.0 - h / 2.0, (int)l, h);

  }
}


class Regions extends GeneElements {  
  // Constructor
  Regions() {
    super();  
  }
  
  
  // Draw the regions
  void Draw() {
    noStroke();    
    fill(0.0, 0.5, 0.0);
    
    super.Draw();
  }
  
  
  // Implement the base class AddElement() method
  void AddElement(String values[], HashMap<String, Integer> columnMap) {    
    int regionStartIndex = columnMap.get("region_start");
    int regionEndIndex = columnMap.get("region_end");
    int noteIndex = columnMap.get("note");
    int nameIndex = columnMap.get("feature_type_type_name");
    
    // Get the start and end positions of this region.
    int regionStart = int(values[regionStartIndex]);
    int regionEnd = int(values[regionEndIndex]);
      
    // Get other information for the region
    String note = values[noteIndex];
    String name = values[nameIndex];
      
    // Add a region
    // In plotit.py min and max are used to sort the regionStart and regionEnd values, perhaps because of transcript direction? 
    elements.add(new Region(min(regionStart, regionEnd), max(regionStart, regionEnd), note, name));
  } 
}


class Region extends GeneElement {  
  // Region information
  String note;
  String name;
  
  
  // Constructor
  Region(int _start, int _end, String _note, String _name) {
    super(_start, _end);
  
    note = _note;
    name = _name;  
  }
 
 
  // Draw the region
  void Draw() {
    float h = 2.0;
//    rect(start, height / 2.0 - h / 2.0 + 20, length, h);
    float x = geography.TranscriptToScreen(start, intronScale);
    
    // XXX: Problem here in that regions may span introns...
    rect(x, height / 2.0 - h / 2.0 + 20, geography.TranscriptToScreen(end, intronScale) - x, h);
  }
}


class Results extends GeneElements {  
  // Constructor
  Results() {
    super(); 
  }
  
  
  // Implement the base class AddElement() method
  void AddElement(String values[], HashMap<String, Integer> columnMap) { 
    int hgvsTranscriptIndex = columnMap.get("hgvstranscript");
    int variantEffectIndex = columnMap.get("variant_effect");
    int tagIndex = columnMap.get("tag");
    
    // HGVS is the Human Gene Variation Society     
    String hgvsTranscript = values[hgvsTranscriptIndex];
      
    // Skip invalid hgvs transcript lines
    if (hgvsTranscript.equals("") || hgvsTranscript.equals("?") || hgvsTranscript.equals("None")) return;
      
    String variantEffect = values[variantEffectIndex];
    String tag = values[tagIndex];
      
    // Add a region
    elements.add(new Result(hgvsTranscript, variantEffect, tag));
  }
}


// XXX: Should probably be called Variant
class Result extends GeneElement {
  // Descriptions
  String hgvsTranscript;
  String variantEffect;
  String tag;
  
  
  // Constructor
  Result(String _hgvsTranscript, String _variantEffect, String _tag) {
    // Don't know the start and end yet, so call the no-argument constructor
    super();
    
    hgvsTranscript = _hgvsTranscript;
    variantEffect = _variantEffect;
    tag = _tag;
  
    SetPosition();
  }
  
  
  // Parse the hgvs transcript string to set the position information
  void SetPosition() {    
    // Extract the position from the hgvs transcript string.
    String s1[] = split(hgvsTranscript, ':');
    String s2[] = split(s1[1], '.');
    
    String s = s2[1];
    
    int i;
    for (i = 0; i < s.length(); i++) {
      if (Character.isLetter(s.charAt(i)) || s.charAt(i) == '_') {
        break;
      }
    }
        
    start = int(s.substring(0, i));
    
    if (s.charAt(i) == '_') {
      // Find the end
      int j;
      for (j = i + 1; j < s.length(); j++) {
        if (Character.isLetter(s.charAt(j))) {
          break;
        }
      }
      
      end = int(s.substring(i + 1, j));
    }
    else {
      end = start;
    }
    
    ComputeLength();
  }
  
  
  // Draw the variant
  void Draw() {     
    float x = (float)(end + start) / 2.0;
    
    
    x = geography.TranscriptToScreen((int)x, intronScale);
    
    
    float glyphScale = 10.0;
    
    float opacity = 0.5;
    
    if (variantEffect.equals("nonsense")) {
      fill(1.0, 0.5, 0.0, opacity / 2.0);
      
      strokeWeight(1);
      stroke(0.75, 0.0, 0.0, opacity);
      
      float y = height / 2.0 + 40;
      
      ellipseMode(CENTER);
      ellipse(x, y, glyphScale / sceneScale, glyphScale);      
    }
    else if (variantEffect.equals("frameshifting indel")) {
      strokeWeight(1);
      stroke(0.75, 0.0, 0.0, opacity);
      
      float y = height / 2.0 + 40;
      float w = glyphScale / 2.0 / sceneScale;
      float h = glyphScale / 2.0;
      
      line(x - w, y - h, x + w, y + h);
      line(x - w, y + h, x + w, y - h);
    }
    else if (variantEffect.equals("missense") && tag.equals("DM")) {      
      stroke(0.75, 0.0, 0.75, opacity);
      strokeWeight(2);
      
      float y = height / 2.0 + 20;
      float w = glyphScale / 2.0 / sceneScale;
      float h = glyphScale / 2.0;
      
      line(x, y - h, x, y + h);
//      line(x - w, y, x + w, y);
    }
  }
}


class Hotspots extends GeneElements {    
  // Constructor
  Hotspots() {
    super();  
  }
  
  
  // Draw the hotspots
  void Draw() {       
    noStroke();
    fill(0.0, 0.0, 1.0);
    
    super.Draw();
  }
    
    
  // Implement the base class AddElement() method
  void AddElement(String values[], HashMap<String, Integer> columnMap) { 
    int pValueIndex = columnMap.get("p_value");
    int permutationsIndex = columnMap.get("N_perm");
    int hotSpotStartIndex = columnMap.get("start");
    int hotSpotEndIndex = columnMap.get("end");
        
    float pValue = float(values[pValueIndex]);
    int permutations = int(values[permutationsIndex]);
    int hotSpotStart = int(values[hotSpotStartIndex]);
    int hotSpotEnd = int(values[hotSpotEndIndex]);
      
    // Add a hotspot
    elements.add(new Hotspot(hotSpotStart, hotSpotEnd, pValue, permutations));
  }
}


class Hotspot extends GeneElement {
  // Descriptions
  float pValue;
  int permutations;
  
  
  // Constructor
  Hotspot(int _start, int _end, float _pValue, int _permutations) {
    super(_start, _end);
    
    pValue = _pValue;
    permutations = _permutations;
  }
  
  
  // Draw the hotspot
  void Draw() {
    float h_height = 40.0;
    float h_factor = 40.0;
    
    float h = 2.0;
    
    float p;
    if (pValue != 0) {
      p = pValue;
      fill(0.0, 0.0, 1.0);
    }
    else {
      float pv = 1.0 / (float)(permutations + 1);
      float y = h_height + h_factor * log(pValue) / log(10);
    }
      
    float y = h_height + h_factor * log(pValue) / log(10);
//    rect(start, height / 2.0 - h / 2.0 - y, length, h);
    float x = geography.TranscriptToScreen(start, intronScale);
    rect(x, height / 2.0 - h / 2.0 - y, length, h);
  }
}
