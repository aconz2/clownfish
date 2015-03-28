
.PHONY: all

files=genes.occur.raw genes.occur genes.occur.raw.rt genes.occur.raw.clamped
counter=clownfish count --genes genes.fa --samples reads.fa --kmer-length 3 --size 1000 


#all: $(files)
all: genes.occur.raw

clean:
	rm $(files)

genes.occur.raw:
	$(counter) > $@

