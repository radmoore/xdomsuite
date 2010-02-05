#!/usr/bin/ruby -w

require 'xdomsuite'

xdom = XDOM.new(ARGV[0], 0.001)

puts xdom.res_coverage
puts xdom.prot_coverage
#xdom.arr_dist.each {|k, v| puts "#{k}, #{v}"}
puts xdom.arr_dist.inspect
puts xdom.collapse.to_s
#xdom.arr_dist.each {|e| puts e}
exit
while(xdom.has_next)

  p = xdom.next_prot
  next if p.has_domains?
  puts p
  exit

end




