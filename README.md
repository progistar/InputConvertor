# InputConvertor

## Quick usage to make a comprehensive proteogenomic sequence database
Build a NeoDB from StringTie, Arriba and variant calls.<br>
```bash
java -Xmx12G -jar makeNeoDB.jar \
--stringtie [StringTie gene annotation file path] \
--fpkm [FPKM threshold for the transcripts predicted by StringTie] \
--genome [Reference genome fasta file path] \
--var [Variant call file path] \
--arriba [Arriba TSV file path] \
--output [Output file path]
```

And then, remove overlapped sequences from filtering sequence databases (it can be multiple fasta files from UniProt and GENCODE).<br>
```bash
java -Xmx12G -jar refineNeoDB.jar \
--reference [A list of reference fasta files (separated by comma). Theses records will be annotated as PE=1] \
--non-reference [A list of non-reference fasta files (e.g., NeoDB.fasta, separated by comma). Theses records will be annotated as PE=2] \
--filter [A list of protein sequences to be filtered in non-reference fasta files (separted by comma).] \
--output [Output file path]
```
