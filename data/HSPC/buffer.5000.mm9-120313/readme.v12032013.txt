All network information can be extracted by running: 

  `$ tar -xzf networks.v12032013.tgz`

Each directory (cell type, n=25) includes a file (genes-regulate-genes.txt)
containing all of the TF-to-TF regulatory interactions mapped within that cell
type.

Regulatory interactions were identified in each cell type by scanning the
proximal DHSs (+/-5kb from the canonical transcriptional start site) of each
transcription factor gene DNaseI footprints corresponding to the recognition
sequences of known TFs.

The file genes-regulate-genes.txt contains two columns of gene names, at column 
4 and column 5 of the text file.

The gene in the first column is bound (i.e., contains a DNaseI footprint in 
proximal regulatory DNA) corresponding to the second column's TF gene product.

The final column (column 6) represents the number of times a motif for the 
regulating gene is found in the regulatory region (promoter) for the regulated 
gene.
