/*=========================================================================
 
 Name:        Geography.pde
 
 Author:      David Borland, The Renaissance Computing Institute (RENCI)
 
 Copyright:   The Renaissance Computing Institute (RENCI)
 
 Description: Contains geography of different sequences for a gene and 
              methods for transforming between coordinate systems.
 
 =========================================================================*/


class Geography {  
  // XXX: There should be four coordinate systems, genomic, relative genome--gene (maybe? screen?), transcript, and coding
  
  int geneOffset;
  int codingOffset;
  
  private class Sequence {      
    // The start and end positions in gene (relative genomic) coordinates
    int start;
    int end;
    
    int intronOffset;
    
    // XXX: Enumerations are not currently supported in Processing...   
    public static final int UTR = 0;
    public static final int EXON = 1;
    public static final int INTRON = 2;
    
    int type;
    
    Sequence(int _start, int _end, int _intronOffset, int _type) {
      start = _start;
      end = _end;
      intronOffset = _intronOffset;
      type = _type; 
    }
  }
  
  ArrayList<Sequence> sequences;
  
  
  Geography() {
    sequences = new ArrayList();  
  }
  
  
  float GetScreenWidth(float intronScale) {
    return GeneToScreen(sequences.get(sequences.size() - 1).end, intronScale);
  }
  
  
  int GenomicToGene(int genomic) {
    return genomic - geneOffset;  
  }
  
  int GeneToGenomic(int gene) {
    return gene + geneOffset; 
  }
  
  
  float GeneToScreen(int gene, float intronScale) {
    // XXX: Could binary search here..
    for (int i = 0; i < sequences.size(); i++) {
      Sequence s = sequences.get(i);
      if (gene >= s.start && gene <= s.end) {
          return gene - s.intronOffset * (1.0 - intronScale);
      }
    }
    
    return -1;  
  }
  
  float TranscriptToScreen(int transcript, float intronScale) {
   /*
      for (int i = 0; i < sequences.size(); i++) {
      Sequence s = sequences.get(i);
      if (s.type != Sequence.INTRON) {
        float screenPos = (float)transcript + (float)s.intronOffset * intronScale;
        if (screenPos >= s.start && screenPos <= s.end) {
          return screenPos;
        }
      }
    }
   
    return -1;  
   */
    return GeneToScreen(TranscriptToGene(transcript), intronScale); 
  }
  
  
  int GeneToTranscript(int gene) {
    // XXX: Could binary search here..
    for (int i = 0; i < sequences.size(); i++) {
      Sequence s = sequences.get(i);
      if (gene >= s.start && gene <= s.end) {
        if (s.type == Sequence.INTRON) {
          // XXX: Should be able to return in N+1 mode as described by Chris Bizon
          return -1; 
        }
        else {
          return gene - s.intronOffset;
        }
      }
    }
    
    return -1;
  }
  
  int TranscriptToGene(int transcript) {
    for (int i = 0; i < sequences.size(); i++) {
      Sequence s = sequences.get(i);
      if (s.type != Sequence.INTRON) {
        int gene = transcript + s.intronOffset;
        if (gene >= s.start && gene <= s.end) {
          return gene;
        }
      }
    }
   
    return -1; 
  }
/*
  int GenomicToTranscript(int genomic) {
    return genomic;
  }
  
  int TranscriptToGenomic(int transcript) {
    for (int i = 0; i < sequences.size(); i++) {
      if (
    }
    
    return transcript;
  }
*/
/*
  double GenomicToScreen(int genomic, double intronFraction) {
    // XXX: Could binary search here..
    for (int i = 0; i < sequences.size(); i++) {
      if (transcript >= sequences[i].start && transcript <= sequences[i].end) {
        if (sequences.type == INTRON) {
          // XXX: Should be able to return in N+1 mode as described by Chris Bizon
          return -1.0; 
        }
        else {
          return transcript + sequences[i].intronSum * intronFraction;
        }
      }
    } 
  }

  double TranscriptToScreen(int transcript, double intronFraction) {
    // XXX: Could binary search here..
    for (int i = 0; i < sequences.size(); i++) {
      if (transcript >= sequences[i].start && transcript <= sequences[i].end) {
        if (sequences.type == INTRON) {
          // XXX: Should be able to return in N+1 mode as described by Chris Bizon
          return -1.0; 
        }
        else {
          return transcript + sequences[i].intronSum * intronFraction;
        }
      }
    } 
  }
*/

  int TranscriptToCoding(int transcript) {
    return max(-1, transcript - codingOffset);
  }
  
  int CodingToTranscript(int coding) {
    return max(-1, coding + codingOffset);
  }  
  
  
  // Read the coordinate data from a geography text file
  void ReadFile(String fileName) {
    // Load the file into an array of strings
    String lines[] = loadStrings(fileName);
    

    // Looks like the coding start and end are always the two values on line 6 in the file.
    // Also looks like tabs are always the delimiter.
    int codingIndex = 5;
    String codingLine[] = split(lines[codingIndex], '\t');

    // Do a simple test until we have a better specification
    if (codingLine.length == 2) {      
      codingOffset = int(codingLine[0]) - 1;
    }
    else {
      println("Coordinates::ReadFile() : Coding start and end not where expected");
      println(codingLine);

      return;
    } 


    // Looks like the transcript start is always 1
/*
    transcriptStart = 1;

    // Looks like the transcript end is the last number on the last line
    String transcriptLine[] = split(lines[lines.length - 1], '\t');
    transcriptEnd = int(transcriptLine[4]);

    transcriptLength = transcriptEnd - transcriptStart + 1;
*/    
    
    // Looks like the transcript starts on the line after the coding start and end
    int geneStartIndex = codingIndex + 1;
    String geneStartLine[] = split(lines[geneStartIndex], '\t'); 
    geneOffset = int(geneStartLine[0]) - 1;
    
    
    // Parse the sequences
    int intronOffset = 0;
    for (int i = geneStartIndex; i < lines.length; i++) {
      // Split the current line by tabs
      String sequenceLine[] = split(lines[i], '\t');
      
      // Get the start and end positions of this sequence in gene coordinates.
      int sequenceStart = int(sequenceLine[0]) - geneOffset;
      int sequenceEnd = int(sequenceLine[1]) - geneOffset;
      
      // Get the type
      String typeString = sequenceLine[2];
      int type;
      if (typeString.contains("UTR")) { 
        type = Sequence.UTR;
      }
      else if (typeString.equals("exon")) {
        type = Sequence.EXON;
      }
      else if (typeString.equals("intron")) {
        type = Sequence.INTRON;
      }
      else {
        println("Coordinates::ReadFile() : Unknown sequence type: " + typeString);
        
        continue;
      }
      
      // Add a sequence 
      sequences.add(new Sequence(sequenceStart, sequenceEnd, intronOffset, type));
      
      
      if (type == Sequence.INTRON) {
        intronOffset += sequenceEnd - sequenceStart + 1; 
      }
    }
  }
}

