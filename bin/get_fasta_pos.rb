#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/util/simpleopt'

$0 = "rbbt #{$previous_commands*" "} #{ File.basename(__FILE__) }" if $previous_commands

options = SOPT.setup <<EOF

get fasta position

$ #{$0} [options] filename.fa contig pos

-h--help Print this help

EOF
if options[:help]
  if defined? rbbt_usage
    rbbt_usage 
  else
    puts SOPT.doc
  end
  exit 0
end

file, contig, pos = ARGV
pos = pos.to_i - 1

Open.open(file) do |f|

  line = f.gets while ! (line && line.include?(contig))
  puts line

  line = f.gets
  puts line[pos-5..pos+5]
  puts line[pos-2..pos+2]
  puts line[pos]
end


