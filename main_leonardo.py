from Bio.Blast import NCBIWWW
from Bio.Blast import Applications
from Bio import SeqIO
help(NCBIWWW.qblast)
record1 = SeqIO.read("GCF_001402945.1_ASM140294v1_genomic.fna", format="fasta")
record1 = SeqIO.read("GCF_001402935.1_ASM140293v1_genomic.fna", format="fasta")
result_handle = NCBIWWW.qblast("blastn", database = "SP1/SP1"  ,sequence=record1.Seq)

def blastn(fas, db):
    print("blasting against the " + db + " database")

    blast_in = (fas)
    blast_out = (tmpdir + "blast_files/" + sample + ".blast")


    blastn_cline = NcbiblastxCommandline(query=blast_in, db=db, evalue=0.001, outfmt=6, out=blast_out, num_threads=12)

    stdout, stderr = blastn_cline()

    blast_err_log = open(tmpdir + "blast_err.txt", "w")
    blast_stdout_log = open(tmpdir + "blast_stdout.txt", "w")

    blast_err_log.write(stderr)
    blast_stdout_log.write(stdout)

    return blast_out
