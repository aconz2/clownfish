
.PHONY: all

files=genes.occur.raw genes.occur genes.occur.raw.rt genes.occur.raw.clamped

all: $(files)

clean:
	rm $(files)

mer_counts.jf: reads.fa
	jellyfish count -m 3 -s 1000 $<

genes.occur.raw: genes.fa mer_counts.jf
	clownfish --genes genes.fa --database mer_counts.jf --raw > $@

genes.occur.raw.rt: genes.fa mer_counts.jf
	clownfish --genes genes.fa --database mer_counts.jf --right-trim 4  --raw > $@

genes.occur.raw.clamped: genes.fa mer_counts.jf
	clownfish --genes genes.fa --database mer_counts.jf --max 1  --raw > $@
	
genes.occur: genes.fa mer_counts.jf
	clownfish --genes genes.fa --database mer_counts.jf --score median > $@

