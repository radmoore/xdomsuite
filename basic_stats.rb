#!/usr/bin/ruby

def mean(values)
  values = string_array2int_array(values)
	values.inject(0) {|sum, x| sum += x} / values.size.to_f
end

# Standard deviation
def stdv(values)
  values = string_array2int_array(values)
	aver = mean(values)
	var = values.inject(0) {|var, x| var += (x - aver)** 2}
	Math.sqrt(var/(values.size-1))
end

# Relative Standard deviation (%RSD)
def rsd(values)
  values = string_array2int_array(values)
  aver = mean(values)
  stdv = stdv(values)
  (stdv*100)/aver
end

# Relative Standard deviation (LOW)
def rel_stdv(values)
  values = string_array2int_array(values)
  aver = mean(values)
  aver/values.length
end

# Z-Score (distance of mean to outliers
# devided by standard deviation - that is,
# number of stdv steps necessary to reach
# outlier)
def z_score(score, values)
  (score.to_i - mean(values))/stdv(values)
end


private

def string_array2int_array(values)
  values.collect{|v| v.to_i}
end

