#!/usr/bin/ruby -w

require 'xdomsuite'

hsap = Proteome.new(ARGV.shift, 0.001)

p1 = hsap.grab
p2 = hsap.grab

puts "P1: ", p1
puts "P2: ", p2

puts "The edit distance from p1 to p2 is: "+ p1.edit_distance(p2).to_s+" operations"

