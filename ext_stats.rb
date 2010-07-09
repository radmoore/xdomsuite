#!/usr/bin/ruby

require 'rubygems'
require 'rsruby'

R = RSRuby.instance

def fisher_test(values, alternative='two.sided')
  R.matrix.autoconvert(RSRuby::NO_CONVERSION)
  m = R.matrix(values,:nrow=>2,:ncol=>2)
  result = R.fisher_test(m, :alternative => alternative)
  pvalues = result["p.value"]
  return pvalues
end

def p_adjust(values, method='fdr')
  result = R.p_adjust(values, :method => method)
  return result
end
