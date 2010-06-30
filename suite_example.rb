#!/usr/bin/ruby -w

# example use of xdomsuite

require 'xdomsuite'

puts "Species, Total proteins, res coverage, prot coverage, total domains, total uniq domains, pfam A domains, pfam B domains, uniq arr (uncollapsed), average dom / protein (uncollapsed), uniq arr (collapsed), average dom / protein (collapsed)"

ARGV.each do |file|
  p = Proteome.new(file, 0.001)
  p.simple_overlap_resolution
  print "#{p.species}, #{p.total_proteins}, #{p.res_coverage}, #{p.prot_coverage}, "
  print "#{p.total_domains.to_s}, #{p.uniq_domains.to_s}, #{p.pfam_A.size.to_s}, #{p.pfam_B.size.to_s},"
  print " #{p.uniq_arrangements.to_s}, #{p.average_domain_no.to_s}, "
  p.collapse!
  print "#{p.uniq_arrangements.to_s}, #{p.average_domain_no.to_s}\n"
end
