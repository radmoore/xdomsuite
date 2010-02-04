#!/usr/bin/ruby -w

require 'xdomsuite'

xdom = XDOM.new(ARGV[0], 0.001)

puts xdom.res_coverage
puts xdom.prot_coverage
exit
while(xdom.has_next)

  p = xdom.next_prot
  next if p.has_domains?
  puts p
  exit

end




