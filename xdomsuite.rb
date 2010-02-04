#!/usr/bin/ruby -w
##!/usr/bin/ruby1.9 -w

# XDOMSUITE
# rv 2.4a

# TODO
# * Use syteny to refine fusion candidates ->
# * include Enumerable (and overwrite each etc)
# * Marshalling?
# * extend to allow hmmoutput
# * if protein length unknown, take stop pos of last domain as length

require 'ftools'

# Add to hash method to array
class Array
	def to_h
   	inject({}) { |m, e| m[e[0]] = e[1]; m }
  end
end



# CLASS: PARSER
class Parser

  # TODO
  # add native hmmer types
  def initialize(filename)
    @filename = filename
    @file_ext = File.extname(@filename)
    @file_bn = File.basename(@filename, @file_ext)
    @comment = '#'
    @hits = Hash.new
  end

  attr_reader :filename, :file_ext, :file_bn
  attr_accessor :comment

  # <length> <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> 
  # <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> <predicted_active_site_residues>
  # * custom added
  def pfamscan(evalue=10, name=true)
    hmmout = File.open(@filename, "r")
    p = nil
    while(line = hmmout.gets)
      next if (/^#{@comment}/.match(line) || /^$/.match(line))
      line.chomp!
      (seq_le, # sequence length, custom field
			 seq_id, # seq  id
       aln_st, # alignment start
       aln_en, # alignment end
       env_st, # envelope start
       env_en, # envelope end
       hmm_ac, # hmm acc
       hmm_na, # hmm name
       dom_ty, # type
       hmm_st, # hmm start
       hmm_en, # hmm end
       hmm_ln, # hmm length
       bit_sc, # bit score
       eva_ht, # e-value
       sig_ht, # significance
       cla_id, # clan
       pre_as) = line.split  # predicted_active_site_residues

      next if (eva_ht.to_f >= evalue)
      name = (name) ? hmm_na : hmm_ac
      @hits[seq_id] = Protein.new(seq_id, seq_le.to_i) unless(@hits.has_key?(seq_id))
      p = @hits[seq_id]
      d = Domain.new(env_st.to_i, env_en.to_i, name, eva_ht.to_f, p.pid)
      p.add_domain(d)

    end
    hmmout.close
    return @hits.values
  end

  private

  # TODO
  def validate(type)
    

  end


end

# CLASS: XDOM
class XDOM

  HEADERre = /^>(\S+)\s+(\d+)$/
  FASTAHEAD = /^>(\S+)\s+.+/
  include Enumerable

  def initialize(filename, evalue, species=String.new)
    @filename = filename
    @species  = (species.empty?) ? get_species(filename) : species
    @evalue   = evalue
    @proteins = Hash.new
    @arrangements = Hash.new
    @domains  = Hash.new
    @total_domains    = 0
    @residue_coverage = 0.00
    @total_dom_residues = 0
    @total_prot_length = 0
    read_xdom(filename)
    @total_proteins   = @proteins.size
    @uniq_domains   = @domains.keys.size
    @current_prot = 0
  end
  attr_accessor :filename, :evalue, :species, :total_proteins, :total_domains, :uniq_domains

  def to_s
    lines = String.new
    @proteins.values.each {|p| lines += p.to_s}
    return lines
  end

	def uniq_arrangements
		@arrangements.keys.size
	end

  def arrangements
    @arrangements.keys
  end

	def get_all_pfam_A
		pfamA = Array.new
		@proteins.values.each {|p| pfamA << p.pfam_A unless p.pfam_A.empty?}
		return pfamA.flatten!
	end

	def get_all_pfam_B
		pfamB = Array.new
		@proteins.values.each {|p| pfamB << p.pfam_B unless p.pfam_B.empty?}
		return pfamB.flatten!
	end

	def has_next
		return (@current_prot == self.total_proteins) ? false : true
	end

  def next_prot
    pid = @proteins.keys[@current_prot]
    @current_prot += 1
    return @proteins[pid]
  end

	def get_collapsed_arrangements
		arr = Array.new
		@proteins.values.each {|p| arr << p if p.collapsed?}
		return arr
	end

  def has_prot?(pid)
    @proteins.has_key?(pid)
  end

  def rewind
    @current_prot = 0
     return
  end

  def domains(did="")
    alldoms = Array.new
    if (did.empty?)
      @domains.values.each {|pid| alldoms.concat(self.get_prot(pid).domains)}
    else
      @domains[did].each{|pid| alldoms.concat(self.get_prot(pid).find_domains(did))}
    end
    return alldoms
  end

	def filter_by_type (type)
		prot_filt = Array.new
		@proteins.values.each {|p| prot_filt << p.type_filter(type)}		
		return prot_filt
	end

	def filter_by_type! (type)
		@proteins.values.each {|p| p.type_filter!(type)}		
	end

	def dom_types
		types = Hash.new
		@proteins.values.each do |p|
			p.domains.each {|d| types[d.type] = (types.has_key?(d.type)) ? types[d.type].succ : 1 }
		end
		return types.sort{|a, b| b[1] <=> a[1]}.to_h
	end


	def average_domain_no
		( (@proteins.values.inject(0){ |sum, p| sum += p.domains.size.to_f}) / @proteins.size.to_f )
	end

  def collapse (repunits=2)
    prot_coll = Array.new
    @proteins.each do |id, prot|
      prot_coll << prot.collapse(repunits)
    end
    return prot_coll
  end

  def collapse! (repunits=2)
    @proteins.values.each {|p| p.collapse!(repunits)}
		update_arrangements()
		self.rewind
  end

  def attach_fasta(fastafile)
    f = File.open(fastafile, "r")
    seq = String.new
    prot = false
    while(line = f.gets)
      line.chomp!
      if (m = FASTAHEAD.match(line))
        if (prot)
          @proteins[prot.pid].sequence = seq
          prot.domains.each{|d| d.sequence = seq[(d.from-1)..(d.to-1)]}
        end
        seq = ""
        prot = @proteins.fetch(m[1], false)
        next
      end
      next unless (prot)
      seq += line
    end
    if (prot)
      @proteins[prot.pid].sequence = seq
      prot.domains.each{|d| d.sequence = seq[(d.from-1)..(d.to-1)]}
    end
    f.close
  end

  # DONE
  def attach_location(gffile)
    prev_line = String.new
    prev_loc = String.new
    gff = read_gff(gffile)

    gff.keys.each do |strand|
      gff[strand].each do |line|

        fields = line.split
        loc     = fields[0]
        from    = fields[3]
        to      = fields[4]
        cline   = fields[8] + ";#{loc};#{from};#{to};#{strand}"
        pid     = fields[8].split(';')[0].split('=')[1]
        p = self.get_prot(pid)
        darr    = p.arr_str

        p.set_location(loc, strand, from, to)

        unless (prev_line.empty?)

          unless (loc == prev_loc)
            prev_line = cline
            prev_loc  = loc
            next
          end

          ppid = prev_line.split(';')[0].split('=')[1]
          # current protein downstream of previous protein
          self.get_prot(ppid).set_downstream(p.pid)
          # previous protein upstream of current protein
          p.set_upstream(ppid)

        end
        prev_line = cline
        prev_loc  = loc
      end
    end
  end

  def get_prot (pid)
    return (@proteins.member?(pid)) ? @proteins[pid] : false
  end

  def has_arr? (arrstr)
    return (@arrangements.keys.member?(arrstr))
  end


  def find_prot_by_arr (arrstr)
      return (@arrangements.member?(arrstr)) ? @arrangements[arrstr] : false
  end

  def find_prot_by_dom (did)
    prots = Array.new
    (@domains[did].collect{|d| d.pid}).uniq.each do |pid|
      prots << self.get_prot(pid)
    end
    return prots
  end

  # TODO
  def find_rsp_cand (xdom2, strict=false)

    candidates = Array.new
    seen = Hash.new

    @arrangements.keys.each do |a1|
      next if (xdom2.has_arr?(a1))
      darr = a1.split(';')
      next if (darr.length == 1)
      1.upto(darr.length-1) do |i|
        cand1 = darr[0...i].join(';') << ';'
        cand2 = darr[i...darr.length].join(';') << ';'
        # decomposed instances of darr may
        # *NOT* occur in self (very stringent)
        next if (strict && ( self.has_arr?(cand1) || self.has_arr?(cand2) ))
        if ( xdom2.has_arr?(cand1) && xdom2.has_arr?(cand2) )
          candidates.push([cand1, cand2, a1])
        end
      end
    end
    
    return candidates

  end

  def res_coverage
    return (@total_dom_residues.to_f/@total_prot_length.to_f)*100
  end

  def prot_coverage
    protwdoms = (@proteins.values.inject(0) {|sum, p| sum.succ if p.has_domains?})
    return ( (100/@total_proteins.to_f) * protwdoms.to_f )
  end


  # TODO
  def annotate_with_context(xdom)

  end


private

	def update_arrangements
		new_arr = Hash.new
		@proteins.each do |pid, p|
			new_arr[p.arr_str] = Array.new unless (new_arr.has_key?(p.arr_str))
			new_arr[p.arr_str].push(pid)
		end
		@arrangements = new_arr
	end

  def get_species (filename)
    file = File.basename(filename, ".xdom")
    file.split('.')[0]
  end

  def read_gff (gffile)
    f = File.open(gffile, "r")
    gff = Hash.new
    gff['+'] = Array.new
    gff['-'] = Array.new
    while(line = f.gets)
      next if (/^##.+/.match(line))
      line.chomp!
      fields = line.split
      next unless (fields[2] == 'protein')
      pid = fields[8].split(';')[0].split('=')[1]
      next unless (self.has_prot?(pid))
      strand  = fields[6]
      next if (strand == '.')
      gff[strand].push(line)
    end
    f.close
    return gff
  end

  def read_xdom (xdomfile)
    f = File.open(xdomfile, "r")
    protein = 0
      while(line = f.gets)
        line.chomp!
        next if (/^$/.match(line) || /^#/.match(line))
        if (m = HEADERre.match(line)) then
          unless (protein == 0)
            @arrangements[protein.arr_str] = Array.new unless (@arrangements.has_key?(protein.arr_str))
            @arrangements[protein.arr_str].push(protein.pid)
          end
          protein = Protein.new(m[1], m[2].to_i, "", @species, "", @filename)
          @total_prot_length += m[2].to_i
          @proteins[m[1]] = protein
          next
        end
        (from, to, did, evalue) = line.split
        domain = Domain.new(from.to_i, to.to_i, did, evalue.to_f, protein.pid)
        @total_dom_residues += (to.to_i-from.to_i)
        @domains[did] = Array.new unless (@domains.member?(did))
        @total_domains += 1 if (@domains[did].push(protein.pid))
        protein.add_domain(domain)
      end
      @arrangements[protein.arr_str] = Array.new unless (@arrangements.has_key?(protein.arr_str))
      @arrangements[protein.arr_str].push(protein.pid)
    f.close
  end

end

# CLASS: PROTEIN
# TODO
# * attr_reader instead of accessor
class Protein

  FASTAHEAD = /^>(\S+)\s+.+/

  def initialize (pid, length = 0, sequence=String.new, species=String.new, comment=String.new, xdom=String.new)
    @pid       = pid
    @length   = length
  	@sequence = sequence
    @species  = species
    @comment  = comment
    @domains  = Array.new
    @arr_str  = String.new
    @current_dom = 0
    @position = Hash.new
		@deleted = nil
		@collapsed = false
    self.set_downstream(false)
    self.set_upstream(false)
    self.set_location(nil, nil, nil, nil)
  end
  attr_accessor :pid, :length, :sequence, :species, :comment, :arr_str
	attr_reader :deleted

  def add_domain (domain)
    @domains.push(domain)
    update_arrangement(domain.did)
  end

  def domains
    return @domains
  end

  def has_domains?
    (@domains.empty?) ? false : true
  end

  def total_domains
    @domains.size
  end

	def uniq_doms
		uniq = Hash.new
		(@domains.each {|d| uniq[d] = (uniq.has_key?(d)) ? uniq[d].succ! : 1 }).size
	end

	def pfam_A
		p = self.type_filter('A')
		return p.domains
	end

	def pfam_B
		p = self.type_filter('B')
		return p.domains
	end

	def type_filter (type)
		p = Protein.new(@pid, @length, @sequence, @species, @comment)
		@domains.each {|d| p.add_domain(d) if (d.type == type)}
		return p
	end

	def type_filter! (type)
		doms_filtered = Array.new
		@domains.each {|d| doms_filtered << d if (d.type == type) }
		@domains = doms_filtered
		return self
	end

  def member? (did)
    @domains.each {|d| return true if (d.did == did)}
    return false
  end

	def collapsed?
		return @collapsed
	end


  def set_location (chromosome, strand, from, to)
    @position['location'] = Hash.new unless (@position.has_key?('location'))
    @position['location']['chomosome'] = chromosome
    @position['location']['strand'] = strand
    @position['location']['from'] = from
    @position['location']['to'] = to
    return
  end

  def set_upstream (pid)
    @position['upstream'] = pid
    return
  end

  def set_downstream (pid)
    @position['downstream'] = pid
    return
  end

  def location
    return @position['location']
  end

  def get_chromosome
    return @position['location']['chomosome']
  end

  def get_strand
    return @position['location']['strand']
  end

  def upstream_arr

  end

  def downstream_arr

  end

	# TODO: return true if overlaps exist
  def overlaps?
    olap = 0
    cdom
    self.domains.each do |d|
      

    end
  end

  # TODO
  # HACK alarm !
  # + sanitize
  # + shorten, check effciency
  # ACON: preference pfamA, maximize evalue
  # ACOV: preference pfamA, maximize coverage
  # BCON: preference pfamB, maximize evalue
  # BCOV: preference pfamB, maximize coverage
  def resolve_overlaps(mode=nil)
    raise "Please select overlap resolution mode: AE=0, AC=1, BE=2, BC=3" if (mode.nil?)
    raise "Invalid resolution mode: AE=0, AC=1, BE=2, BC=3" unless ((0 <= mode) && (mode <= 3))

    pos = 0
    pdom = nil
    deleted = Array.new
    until (pos == self.total_domains)
      cdom = self.domains[pos]
      if ( pdom.nil? || (pdom == cdom) || (cdom.from > pdom.to) )
        pdom = cdom
        pos += 1
        next
      end
      unless (cdom.type == pdom.type)
        if (((mode == 0)||(mode == 1)))
          if (cdom.type == 'B')
            deleted.push(cdom)
            self.domains.delete(cdom)
          else
            deleted.push(pdom)
            self.domains.delete(pdom)
          end
         elsif ( (mode == 2)||(mode == 3) )
          if (cdom.type == 'A')
            deleted.push(cdom)
            self.domains.delete(cdom)
          else
            deleted.push(pdom)
            self.domains.delete(pdom)
          end
        end
        pos = 0
        pdom = nil
        next
      end
			# in mode 1, 3: if length is equal, no conflict resolution possible
			# than resolve by evalue anyways
      if (mode%2 == 0 || (cdom.length == pdom.length))
        if (cdom.evalue >= pdom.evalue)
          deleted.push(cdom)
          self.domains.delete(cdom)
        else
          deleted.push(pdom)
          self.domains.delete(pdom)
        end
      else
        if (cdom.length < pdom.length)
          deleted.push(cdom)
          self.domains.delete(cdom)
        else
          deleted.push(pdom)
          self.domains.delete(pdom)
        end
      end
      pos = 0
      pdom = nil
    end
		# if there are interdomainic regions...
		unless (self.get_interdoms.keys.length == 0)
    	pdom = nil
    	self.domains.each_index do |pos|
      	cdom = self.domains[pos]
      	unless (pdom.nil?)
					# try to fit some of the removed domains
					# into the final arrangement
      		deleted.each do |remdom|
	       		next unless ((remdom.from > pdom.to) && (remdom.to < cdom.from))
	        	self.domains.insert(pos, remdom)
	        	deleted.delete(remdom)
	      	end
    		end
				pdom = cdom
			end
			# TODO
			# this should be tried more than once
			# attempt to drop a removed domain into the last position
   		deleted.each do |remdom|
   			next unless ((remdom.from > pdom.to) && (remdom.to < @length))
	    	self.domains.push(remdom)
	      deleted.delete(remdom)
			end
			@deleted = deleted
		end
		return
  end

  def find_domains (did)
    doms = Array.new
    return doms unless self.member?(did)
    @domains.each do |d|
      doms.push(d) if (d.did == id)
    end
    return doms
  end

  def to_s
    head = (@length != 0) ? ">"+@pid+"\t"+@length.to_s+"\n" : ">"+@pid+"\t"+"\t"+"\n"
    doms = String.new
		# NOTE: if protein was created from hmmout, it is possible that the domains
		# (in partcular for e10) and not _quite_ properly sorted. This methods sorts
		# in place, and will repeat the sort everytime its called
		# (***UNEFFCIENT***), to ensure that xdoms are *always* properly sorted
		self.domains.sort! {|x, y| x.from <=> y.from}
    self.domains.each do |d|
      doms += "#{d.from.to_s}\t#{d.to.to_s}\t#{d.did}\t#{d.evalue.to_s}"
      doms += (d.comment.empty?) ? "\n" : "\t;#{d.comment}\n"
    end
    head+doms
  end

  def next_dom
    dom = @domains[@current_dom]
    @current_dom += 1
    return dom
  end

  def rewind
    @current_dom = 0
    return
  end

  def collapse (repunit=2)
    repno = 1
    prev_dom = -1
    rep_doms = Array.new
    doms = Array.new
    self.domains.each do |d|
      unless (prev_dom == -1)
        if (prev_dom.did == d.did)
          repno += 1
          rep_doms.push(prev_dom)
        else
          if (repno >= repunit)
						@collapsed = true
            c_dom = Domain.new(rep_doms[0].from, prev_dom.to, prev_dom.did, -1, self.pid, " #{repno} collapsed repeat")
            doms.push(c_dom)
            rep_doms.clear
          elsif (repno > 1)
            doms += rep_doms.push(prev_dom)
            rep_doms.clear
          else
            doms.push(prev_dom)
          end
          repno = 1
        end
      end
      prev_dom = d
    end
    if (repno >= repunit)
			@collapsed = true
      c_dom = Domain.new(rep_doms[0].from, prev_dom.to, prev_dom.did, -1, self.pid, "#{repno} collapsed repeats")
      doms.push(c_dom)
    elsif (repno > 1)
      doms += rep_doms.push(prev_dom)
    else
      doms.push(prev_dom)
    end
    return doms
  end

  # collapse in place
  def collapse! (repunit=2)
    self.arr_str = String.new
    doms = self.collapse(repunit)
    doms.each {|d| self.arr_str += "#{d.did};" }
    self.comment += "Repeats collapsed if found"
    @domains = doms
  end

  # return procentage if res = false
  # otherwise number of residues covered
  # by domains
  def coverage(res = false)
    domres = 0
    @domains.each {|d| domres += (d.to - d.from)}
    return (res) ? domres : round( ((domres.to_f/self.length.to_f)*100).to_f )
  end


  # TODO
  def edit_distance (protein, collapse=false)
    arr1 = self.arr_str.split(';').sort
    arr2 = protein.arr_str.split(';').sort
    if (collapse)
      arr1.uniq!
      arr2.uniq!
    end
    operations = (arr1.length > arr2.length) ? (arr1.length-arr2.length) : (arr2.length-arr1.length)
    # if add operations are sufficient to reach arrangement 2,
    # than what arr1 and arr2 have in common makes up for the
    # difference between both arrangements
    return operations if ((arr1 & arr2).length == operations)
    for i in 0..(arr1.length) do
      break if (i > arr2.length)
      operations += 2 unless (arr1[i] == arr2[i])
    end
    return operations
  end

  def jaccard_dist (protein, collapse=true)
    arrangement1 = self.domains.collect {|d| d.did}
    arrangement2 = protein.domains.collect {|d| d.did}
    intersect = Array.new
    if (collapse)
      intersect = arrangement1 & arrangement2
      union = (arrangement1 + arrangement2).uniq
    else
      intersect = true_intersection(arrangement1, arrangement2)
      union = (arrangement1 + arrangement2)
    end
    jaccard_coef = (intersect.size.to_f/union.size.to_f)
    return round(1-jaccard_coef)
  end

  def get_interdoms (min=20)
    return nil if (self.total_domains == 0)
    interdoms = Hash.new
    prev_pos  = 0
    dcount    = 0 # starts with pos _before_ first dom
    @domains.each do |d|
      if (prev_pos == 0)
        interdoms[1] = (d.from-1) unless ((d.from-1) < min)
      else
        interdoms[prev_pos] = (d.from-1) unless ((d.from-1) - prev_pos < min )
      end
      prev_pos = d.to+1
      dcount += 1
    end
    if (dcount >= self.total_domains) # when n = total_doms, we are at last interdom
      interdoms[prev_pos] = @length unless ((@length - prev_pos < min))
    end
    return interdoms
  end

  # TODO
  # * wrap fasta output?
  def get_interdom_seq (min=20, fasta=false)
    return false if (@sequence.empty?)
    interdom_seqs = Array.new
    fastaseq = String.new if fasta
    interdoms = self.get_interdoms(min.to_i)
    interdom_no = 0
    interdoms.keys.sort.each do |start|
      interdom_no += 1
      stop = interdoms[start]
      seq = @sequence[(start-1)..(stop-1)]
      interdom_seqs.push(seq)
      if fasta
        fastaseq += ">#{interdom_no}.#{self.pid}\t#{start}-#{stop}\n"
        fastaseq += seq+"\n"
      end
    end
    return (fasta) ? fastaseq : interdom_seqs
  end

  private

  # returns the number of domain overlaps

  def update_arrangement (did)
    @arr_str += "#{did};"
  end

  # do not remove duplicates
  def true_intersection (arr1, arr2)
    arr1.sort!
    arr2.sort!
    intersect = Array.new
    equal = 0
    for i in (0..arr1.length-1)
      for j in (equal..arr2.length-1)
        if arr1[i] == arr2[j]
          equal = j + 1
          intersect[intersect.length] = arr1[i]
          break
        end
      end
    end
    return intersect
  end

  def round(float, pos=2)
    exponent = pos + 2
    float = float*(10**exponent)
    float = float.round
    float = float / (10.0**exponent)
  end

end

# CLASS: DOMAIN
class Domain

  ATYPE = /PF.+/
  BTYPE_ACC = /PB\d+/
# Pfam-B
  BTYPE_ID = /^Pfam-B_.+/
  CTYPE = /^COILS.+/
  DTYPE = /^DIS.+/

  def initialize (from, to, did, evalue, pid, comment=String.new, acc=String.new, go=String.new, sequence=String.new)
    @from = from
    @to = to
    @did = did
    @acc = acc
    @go = go
    @evalue = evalue
    @sequence = sequence
    #@protein = protein
    @pid = pid
    @type = type_check
    @comment = comment
    @length = @to - @from
  end
  attr_accessor :from, :to, :did, :acc, :go, :evalue, :sequence, :type, :pid, :comment, :length

  def fasta_seq
    return (sequence.empty?) ?
      false :
      ">#{self.pid}.#{self.did}\t#{self.from}-#{self.to}\n#{self.sequence}"
  end

  def to_s
    self.from.to_s+"\t"+self.to.to_s+"\t"+self.did+"\t"+self.evalue.to_s
  end

  private

  def type_check
    return 'A' if (ATYPE.match(self.did))
    return 'B' if (BTYPE_ID.match(self.did))
    return 'B' if (BTYPE_ACC.match(self.did))
    return 'C' if (CTYPE.match(self.did))
    return 'D' if (DTYPE.match(self.did))
    return 'A'
  end


end

# CLASS: INDICATOR
class Indicator

  def initialize(indtype = "dots")
    @stop_indicator = false
    @indtype = indtype
  end

  def start
    if (@indtype == "lines")
      loop {
        STDERR.print("\r|")
        STDERR.print("\r/")
        STDERR.print("\r-")
        STDERR.print("\r\\")
        STDERR.print("\r|")
        STDERR.print("\r/")
        STDERR.print("\r-")
        STDERR.print("\r\\")
        break if @stop_indicator
      }
      elsif (@indtype == "dots")
        loop {
          STDERR.puts(".")
          break if @stop_indicator
        }
     end
  end

  def stop
    @stop_indicator = true
    STDERR.puts
  end




end


