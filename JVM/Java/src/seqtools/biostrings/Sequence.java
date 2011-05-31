package seqtools.biostrings;

interface Sequence<T> {
  
  public Iterable<T> getSequence();
  public Sequence<T> view(long start, long end);  
  public long length();
  
}
