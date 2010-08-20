#!/usr/bin/ruby -w

require 'xdomsuite'

hsap = Proteome.new(ARGV.shift, 0.001)

prot = hsap.grab

puts "Searching with #{prot.arrstr}"
test = hsap.collect_by_similarity(1)
test = hsap.find_by_similarity(prot, 2, 5)

puts "Found #{test.length} results:"
test.each{|p| puts p}
