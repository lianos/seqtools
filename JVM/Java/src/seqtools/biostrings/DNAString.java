package seqtools.biostrings;

class DNAString extends Sequence {
  
  public DNAString(String name, byte[] seq) {
    super(name, seq);
  }

  public DNAString(byte[] seq) {
    super("", seq);
  }

}