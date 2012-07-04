```ruby
 
 require 'xdomsuite'

 # create a protein instance using a match requirement of 25% and an evalue
 # threshold of 0.001
 hsap = Proteome.new('hsap.NCBI36.51.hmmout', 0.001, 25); 

 # get proportion of residues found within domain regions (returns 0.63)
 puts hsap.res_coverage

 # get the number of proteins with at least one domain
 puts hsap.prot_coverage

 doms = Hash.new([])
 architectures = Hash.new(0)
 while(hsap.has_next?)            # iterate over all proteins

   p = hsap.next_prot             # get next protein
   next if p.domains.count == 0   # ignore protein if it contains no domains
   architectures[p.arrstr] += 1   # p.arrstr returns a string representation 
                                  # of the domain architecture of p
   p.collapse!                    # collapse! is an in-place routine to 
                                  # collapse consecutive stretches of the 
                                  # same domain
   # create a map of domain accessions to proteins
   p.domains.each{|d| doms[d.acc] << p}
 end 
 ```