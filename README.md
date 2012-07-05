#### Introduction
The xdomsuite is a collection of classes and methods for dealing with domain data. It consists of three 
main classes:

* Proteome
* Protein
* Domain

A Proteome is a collection of protein objects, and contains methods for searching, sorting, counting, 
filtering, comparing, … arrangements in a proteome. It also contains methods for comparing two proteomes 
with each other (set operations). A Protein is a collection of instances of the Domain class, and 
represents a domain arrangement. Similar to the Proteome class contains various methods for comparison etc, 
as well as methods for attaching annotation, genome location etc.

##### Installation
The xdomsuite is available as a gem. To install, add the source and run gem install
```bash
 $ gem sources --add http://iebservices.uni-muenster.de/radmoore/gems/
 $ sudo gem install xdomsuite
```
Depending on your system settings, you may have to use require ‘rubygems’ before you can require xdomsuite.

##### Recognized input formats

* hmmscan
* pfam_scan

While the suite is currently optimized for Pfam (pfam.sanger.ac.uk/), it can be used in conjuntion with 
any domain database. Proteomes can be instantiated from hmmscan or using output of the pfam_scan wrapper. 
For more information on the pfam_scan utility, consult

ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/README

for more information on HMMER, consult

http://hmmer.janelia.org

To this point, HMMER 3.0 (March 2010) and pfam_scan >= v 1.27 output have been tested.

###### xdom format
We recommend using the xdom format. It is the generic IO format used by the xdomsuite and is a easy 
to read, fasta-like representation of a protein with its constituent domains. For example, the 
human protein ENSP00000376776 is represented in xdom format as:

<pre>
 >ENSP00000376776  617
 57  171 DOMON 2.0e-25
 213 341 Cu2_monooxygen  7.5e-43
 360 521 Cu2_monoox_C  2.3e-52
</pre>

ENSP00000376776 consists of three domains (DOMON, Cu2_monooxygen, Cu2_monoox_C) and has a sequence 
length of 617. In the xdom format, domains are sorted in N- to C terminal order. Each domain line 
contains the start and stop position, some form of ID or accession number and an Evalue. This is also 
the format that a call to Proteome.to_s / Protein.to_s will return.

#### Example

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
 
#### Misc

###### Authors
* Andrew D. Moore (radmoore@uni-muenster.de)
* Lothar Wissler (l.wissler@uni-muenster.de)

###### License
Distributes under the same terms as Ruby

###### Warranty
This software is provided “as is” and without any express or implied warranties, including, 
without limitation, the implied warranties of merchantibility and fitness for a particular 
purpose.

###### How do I Reference?
If your journal supports it, please reference either this page (http://github.com/radmoore/xdomsuite) or http://iebservices.uni-muenster.de/radmoore/xdomsuite