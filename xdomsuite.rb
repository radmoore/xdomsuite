#!/usr/bin/ruby -w

#==xdomsuite - A suite of simple procedures for working with domain data
# 
# The xdomsuite contains 2 primary classes, and two 'utility' classes. The way I use these classes, is that I create either a xdom object
# using XDOM or a collection of Protein objects using Parser. That is, I do not use the classes Protein or Domain directly (hence utility).
# 
# But one could, if one wanted.
# 
#===Primary classes
#
#
# * Parser: Parse hmmer output to proteins and domains
# * XDOM: Parse xdom files to proteome-like objects
# 
#===Utility classes
#
# * Protein: Create protein objects
# * Domain: Create domain objects
# 
#==== Usage
#
#===Open Issues /  TODOs
#
# * Use syteny to refine fusion candidates
# * include Enumerable (and overwrite each etc)
# * Marshalling?
# * extend to allow hmmoutput
# * if protein length unknown, take stop pos of last domain as length
#
#===Credits
#Andrew D. Moore (radmoore@uni-muenster.de)
#
#===Version
# rv 2.4a

require 'ftools'

# Add to hash method to array
class Array
	def to_h
   	inject({}) { |m, e| m[e[0]] = e[1]; m }
  end
#  def to_h
#    h = Hash[*self.flatten]
#  end
end



#=== Parser
#Class for creating protein and domain objects directly from hmmer output. 
#Once called, Parser returns an Array of protein objects
#
#==== Usage
# 
# require 'xdomsuite'
# 
# parser = Parser.new('/path/to/file') 
# proteins = parser.pfamscan(0.0001)
# proteins.each {|p| puts p.pid }
# 
#====Open Issues /  TODOs
# 
# * Currently only supports custom pfamscan output 
# * Check file type validity
class Parser

  # Returns an instance of the Parser class.
  def initialize(filename)
    @filename = filename
    @file_ext = File.extname(@filename)
    @file_bn = File.basename(@filename, @file_ext)
    @comment = '#'
    @proteins = Hash.new
  end

  attr_reader :filename, :file_ext, :file_bn
  attr_accessor :comment

  # Calls the pfamscan parser method (new version), and returns an Array of protein objects. Domain must have
  # a score better than +evalue+ before they are included. If +name+ is +true+, then the domain name is used.
  # Otherwise, the domain accession number is used
  # 
  # The fields in the pfam_scan output file are
  # 
  # * length (custom filed to maintain protein length)
  # * seq id
  # * alignment start
  # * alignment end
  # * envelope start
  # * envelope end
  # * hmm acc
  # * hmm name
  # * type 
  # * hmm start
  # * hmm end
  # * hmm length
  # * bit score
  # * E-value
  # * significance
  # * Clan
  # * Predicted_active_site_residues
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
      did = (name) ? hmm_na : hmm_ac
      did = did.split('.')[0] if (/.+\.\d+/.match(did))
      @proteins[seq_id] = Protein.new(seq_id, seq_le.to_i) unless(@proteins.has_key?(seq_id))
      p = @proteins[seq_id]
      d = Domain.new(env_st.to_i, env_en.to_i, did, eva_ht.to_f, p.pid, cla_id)
      p.add_domain(d)
    end
    hmmout.close
    return @proteins.values
  end



  # TODO: complete
  def hmmscan(evalue=10, name=true, custom=true)
      
      hmmout = File.open(@filename, "r")
      p = d = nil
      length = 0
      
      next if (/^#{@comment}/.match(line) || /^$/.match(line))

      while(line = hmmout.gets)
        line.chomp!
        fields = line.split
        length = fields.shift if (/\d+/.match(fields[0]))
        next if (fields[12].to_f >= evalue)
        name = (name) ? fields[6] : fields[5]
        @proteins[fields[0]] = Protein.new(fields[0], length) unless (@proteins.has_key?(fields[0]))
        p = @proteins[fields[0]]
      end

      hmmout.close
        
  end

  private

  # Check file type for validity
  def validate(type)
  end

end

#=== XDOM
#Class for creating a proteome-like object from a xdom file. Provides methods for iteration and filtering, sequence association
#collapsing repeats, retreiving statistics etc 
#
#==== Usage
# 
# require 'xdomsuite'
#
# xdom = XDOM.new(ARGV[0], 0.001)
# puts xdom.res_coverage
# puts xdom.prot_coverage
# 
# while(xdom.has_next)
#  p = xdom.next_prot
#  next unless p.has_domains?
#  p.domains.each {|d| puts d.did}
# end
#
#
#====Open Issues /  TODOs
# 
# * Currently only supports custom pfamscan output 
# * Check file type validity
class XDOM

  HEADERre = /^>(\S+)\s+(\d+)$/
  FASTAHEAD = /^>(\S+)\s+.+/
  include Enumerable


  @@clans = false
  @@names = true

  # Returns a new xdom object containing proteins and domains
  #
  # * +filename+ : path to the xdomfile
  # * +evalue+   : annotation threshold
  # * +species+  : name of the proteome (OPTIONAL)
  def initialize(filename, evalue, species = "")
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

  attr_reader :filename, :evalue, :species, :total_proteins, :total_domains, :uniq_domains

  def self.clans
    @@clans
  end

  # TODO clean this up
  def self.clans=(var)
    unless (var.kind_of?(TrueClass) or var.kind_of?(FalseClass))
      raise "Class variable 'clans' must of type boolean"
    end
    if var
      begin
        require 'pfam_translator'
      rescue LoadError
        STDERR.puts "ERROR: pfam_translator.rb required to map to clans, setting clans to false"
        var = false
      end
    end
    @@clans = var
  end

  def self.names
    @@names
  end

  # method to set class varable names, that indicates
  # whether the xdom file read to create a new xdom
  # object uses names or ids. This is relevant for acc || id to
  # clan mapping
  def self.names=(var)
    unless (var.kind_of?(TrueClass) or var.kind_of?(FalseClass))
      raise "Class variable 'names' must of type boolean"
    end
    @@names = var
  end

  # returns a random protein object
  # or nil
  def grab
    return nil if @proteins.keys.empty?
    return @proteins[@proteins.keys[rand(@proteins.size)]]
  end

  # Returns string representation of self
  def to_s
    lines = String.new
    @proteins.values.each {|p| lines += p.to_s}
    return lines
  end

  # Returns a list of unique arrangements as an array of strings
	def uniq_arrangements
		@arrangements.keys.size
	end

  # Returns a list of all arrangements as an array of strings
  def arrangements
    @arrangements.keys
  end

  # Returns an array of domain objects
  # 
  # * If _did_ is specified, an array of instances of _did_ is returned, or an empty array if non were found
  # * If _type_ is specified, an array of instances of Domains of _type_ is returned, or an empty array if non were found
  # * If _did_ AND _type_ are specified, method behaves as if only _did_ was specified
  def domains(did = nil, type=nil)
    doms = Array.new
    if (did.nil?)
      if type.nil?
        @proteins.values.each {|p| doms << p.domains}
      elsif (type == 'A')
  		  @proteins.values.each {|p| doms << p.pfam_A unless p.pfam_A.empty?}
      elsif (type == 'B')
  		  @proteins.values.each {|p| doms << p.pfam_B unless p.pfam_B.empty?}
      end
    else
      @proteins.values.each {|p| doms << p.find_doms(did)}
    end 
    return doms.flatten
  end


  # Returns an array of domain identifiers representing the unique set of domains present in self 
  def get_all_uniq_doms
    @domains.keys
  end

  # Intended for iteration over proteins in self. Returns _true_ if self
  # still has a protein, otherwise _false_
	def has_next?
		return (@current_prot == self.total_proteins) ? false : true
	end

  # Returns the next Protein, or _nil_ if no protein is left in the current
  # iteration. Proteins are hashed for fast access, hence the order in which the
  # Protein objects are returned can vary. Uses an internal instance variable to
  # maintain the current position. The position can be reset using rewind
  def next_prot
    pid = @proteins.keys[@current_prot]
    @current_prot += 1
    return @proteins[pid]
  end

  # Resets the position in the current interation to 0 (see next_prot and has_next?)
  def rewind
    @current_prot = 0
    return
  end

  # Returns an array of Proteins with collapsed repeats, or an empty array if non exist
  # ( see Protein.collapse )
	def get_collapsed_arrangements
		arr = Array.new
		@proteins.values.each {|p| arr << p if p.collapsed?}
		return arr
	end

  # Return _true_ if _self_ contains a Protein with the id _pid_, otherwise _false
  def has_prot?(pid)
    @proteins.has_key?(pid)
  end

  # Returns a Hash containing the frequency of proteins with n domains, where _res_ represents the limit.
  # If _collapsed_ is _true_, Proteins are collapsed prior to counting (repeats are not counted).
  #
  # If XDOM contains four Proteins with the arrangements
  #
  #  P1 = [a, b, c, d]
  #  P2 = [a, b]
  #  P3 = [a]
  #  P4 = [a]
  #
  # then arr_dist returns a hash of the following form
  #  { 1 => 2, 2 => 1, 4 => 1 }
  # and if _res_ = 2
  #  { 1 => 2, 2 => 1, '>2' => 1 } 
  #
  # TODO: 
  # * Return sorted hash
  def arr_dist(collapse = true, res=10)
    dist = Hash.new
    @proteins.values.each do |p|
      p = p.collapse if (collapse)
      next unless p.has_domains?
      d_no = p.domains.size
      if (d_no <= res)
        dist[d_no] = (dist.has_key?(d_no)) ? dist[d_no].succ : 1
      else
        k = ">#{res}"
        dist[k] = (dist.has_key?(k)) ? dist[k].succ : 1
      end
    end
    return dist
  end

  # TODO: this does not work - seems to modify self!!!
  # Returns a new xdom where all domain that are not of type _type_ are removed
  # See Protein.type_filter
	def filter_by_type(type)
    newxdom = self.dup
    while(newxdom.has_next?)
      p = newxdom.next_prot
      p.type_filter!(type)
    end
    newxdom.rewind
    return newxdom
	end

  # In place version of filter_by_type 
	def filter_by_type!(type)
		@proteins.values.each {|p| p.type_filter!(type)}		
    update_arrangements()
    return
	end

  # Returns a Hash with the frequency counts for all known _types_ of domains
  # present in _self_
  #
  # Example return:
  #  { 'A' => 4323, 'B' => 32123 }
  #
  # Issues
  # * Not sure if this method is really useful
	def dom_types_freq
		types = Hash.new
		@proteins.values.each do |p|
			p.domains.each {|d| types[d.type] = (types.has_key?(d.type)) ? types[d.type].succ : 1 }
		end
    puts "type before sorting: "+types.inspect
		return types.sort{|a, b| b[1] <=> a[1]}.to_h
	end

  # Returns the average number of domains present in the proteins of self. If _num_ is _true_, returns a float. If _num_ is _false_, returns a formatted string
	def average_domain_no(num = false)
		dno = ( (@proteins.values.inject(0){ |sum, p| sum += p.domains.size.to_f}) / @proteins.size.to_f )
   	return num ? dno : sprintf('%.2f', dno)
	end

  # Returns an Array of Proteins present in _self_, after collapsing repeats. 
  # * repunits : An integer. If a domain is repeated more than _repunits_ times in a row in a protein, it is collapsed to a single instance
  #
  # See Protein.collapse
  def collapse(repunits=2)
    prot_coll = Array.new
    @proteins.each do |id, prot|
      prot_coll << prot.collapse(repunits)
    end
    return prot_coll
  end

  # In place version of collapse
  def collapse!(repunits=2)
    @proteins.values.each {|p| p.collapse!(repunits)}
		update_arrangements()
		self.rewind
  end

  # Provides an interface to attach a FASTA file to the Proteins in _self_
  #
  # Requires that the IDs of the Proteins in _self_ correspond to the IDs in the FASTA file. IDs are matched 
  # using the REGEX described in FASTAHEAD
  # 
  # If a given ID in FASTA can be found in _self_, the sequence in FASTA is set to the corresponding protein using
  # the _sequence_ attribute and, for each domain present in protein, the domain sequence is set using domains _sequence_
  # attribute
  #
  # TODO
  # * return number of successfully attached sequences
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

  # Returns the protein with the ID _pid_, or nil if such a protein does not exist
  def get_prot(pid)
    return (@proteins.member?(pid)) ? @proteins[pid] : nil
  end

  # Returns _true_ if _self_ contains a protein with the arrangement _arrstr_
  # 
  # The _arrstr_ is a string representation of the domains in a protein, where _sep_
  # signifies the character used to seperate the domain IDs.
  # 
  # Example:
  #  "domainA;domainB;domainC", where sep = ';' (default)
  def has_arr?(arrstr, sep=';')
    arrstr.gsub("#{sep}", ';') if (sep != ';')
    return (@arrangements.keys.member?(arrstr))
  end

  # Returns an array of Proteins with the arrangement _arrstr_, or _nil_ if none are found
  #
  # * arrstr : A string representation of the protein domain arrangement
  # * sep : The character that used to seperate the domains in the string
  #
  # See has_arr?
  # 
  def find_prot_by_arr(arrstr, sep=';')
    arrstr.gsub("#{sep}", ';') if (sep != ';')
    pids = (@arrangements.member?(arrstr)) ? @arrangements[arrstr] : nil
    return pids if pids.nil?
    return pids.collect{|pid| @proteins[pid]}
  end

  # Returns an array of proteins that contain at least one domain with the ID _did, or _nil_ if none are found
  def find_arr_by_dom(did)
    prots = Array.new
    return nil unless @domains.has_key?(did)
    (@domains[did].collect{|d| d}).uniq.each do |pid|
      prots << self.get_prot(pid)
    end
    return prots
  end

  # Returns an array of Proteins that contain _did_. If multiple proteins exist that contain a domain with the ID _did_, then only one 
  # (representative) Protein is returned. If _collapse_ is true, proteins are collapsed before determining their uniqueness.
  # 
  # Returns _nil_ if no protein with domain _did_ is present in _self_
  def find_uniq_arr_by_dom(did, collapse=true)
    uniq = Hash.new
    prots = self.find_arr_by_dom(did)
    return nil if prots.nil?
    prots.each do |p|
      p = p.collapse if collapse
      next if uniq.has_key?(p.arr_str)
      uniq[p.arr_str] = p
    end
    return uniq.values
  end

  # Returns the procentage of amino acid residues that fall into an annotated region (domain). 
  # If _num_ is _true_ , a float is returned. Otherwise a formatted string is returned.
  #
  # TODO:
  # * check number return
  def res_coverage(num = false)
    cov = (@total_dom_residues.to_f/@total_prot_length.to_f)*100
    return num ? cov : sprintf('%.2f', cov)
  end

  # Returns the procentage of proteins that contain at least one domain. 
  # If _num_ is _true_ , a float is returned. Otherwise a formatted string is returned.
  #
  # TODO:
  # * check number return
  def prot_coverage(num = false)
    protwdoms = (@proteins.values.inject(0) {|sum, p| sum.succ if p.has_domains?})
    cov = ( (100/@total_proteins.to_f) * protwdoms.to_f )
#    return num ? cov : sprintf('%2f', cov)
  end

  def resolve_overlaps(mode)
    @proteins.values.each {|p| p.resolve_overlaps(mode)}
    return nil
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
        clan = nil
        if @@clans
          if (PfamTranslator::ACC2CLAN.has_key?(did))
            clan = PfamTranslator::ACC2CLAN[did]
          elsif (PfamTranslator::NAME2CLAN.has_key?(did))
            clan = PfamTranslator::NAME2CLAN[did]
          end
        end
        domain = Domain.new(from.to_i, to.to_i, did, evalue.to_f, protein.pid, clan)
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
    @sep = ';'
    @clans = false
    @pfamb = false
  end
  attr_accessor :length, :sequence, :species, :comment, :sep, :clans, :pfamb
	attr_reader :deleted, :pid, :arr_str

  def add_domain (domain)
    @domains.push(domain)
    update_arrstr()
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

  # Wrapper for type_filter. Returns an array of Domain objects of type PfamA, or an empty array if no PfamA domains are found 
	def pfam_A
		p = self.type_filter('A')
		return p.domains
	end

  # Wrapper for type_filter. Returns an array of Domain objects of type PfamB, or an empty array if no PfamB domains are found
	def pfam_B
		p = self.type_filter('B')
		return p.domains
	end

  # Returns a copy of _self_, containing only Domains of _type_. If no Domain objects of _type_ are present in _self_,
  # then self.domains will return an empty array
	def type_filter(type)
		p = Protein.new(@pid, @length, @sequence, @species, @comment)
		@domains.each {|d| p.add_domain(d) if (d.type == type)}
		p.update_arrstr
		return p
	end

  # In place variant of self.type_filter
	def type_filter!(type)
		doms_filtered = Array.new
		@domains.each {|d| doms_filtered << d if (d.type == type) }
		@domains = doms_filtered
		self.update_arrstr
		return self
	end

  def member?(did)
    @domains.each {|d| return true if (d.did == did)}
    return false
  end

	def collapsed?
		return @collapsed
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
		return self
  end

  def find_doms(did)
    doms = Array.new
    return doms unless self.member?(did)
    @domains.each do |d|
      doms.push(d) if (d.did == id)
    end
    return doms
  end

  def get_cooc_doms(did)
    a = Array.new
    self.domains.each {|d| next if (did == d.did); a << d}
    return a
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
      did = d.did
      if @clans
        did = (d.cid.nil?) ? d.did : d.cid
      end
      doms += "#{d.from.to_s}\t#{d.to.to_s}\t#{did}\t#{d.evalue.to_s}"
      doms += (d.comment.empty?) ? "\n" : "\t;#{d.comment}\n"
    end
    head+doms
  end

  def to_arr_s(sep = ';')
    ( self.domains.collect { |d|
        did = d.did
        did = d.cid if (clans && (not d.cid.nil?))
        did
      }
    ).join(sep)
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
    p = Protein.new(@pid, @sequence, @species, @comment)
    self.domains.each do |d|
      unless (prev_dom == -1)
        if (prev_dom.did == d.did)
          repno += 1
          rep_doms.push(prev_dom)
        else
          if (repno >= repunit)
						@collapsed = true
            c_dom = Domain.new(rep_doms[0].from, prev_dom.to, prev_dom.did, -1, self.pid, prev_dom.cid, " #{repno} collapsed repeats")
            p.add_domain(c_dom)
            doms.push(c_dom)
            rep_doms.clear
          elsif (repno > 1)
            doms += rep_doms.push(prev_dom)
            rep_doms.push(prev_dom).each {|dom| p.add_domain(dom)}
            rep_doms.clear
          else
            p.add_domain(prev_dom)
            doms.push(prev_dom)
          end
          repno = 1
        end
      end
      prev_dom = d
    end
    if (repno >= repunit)
			@collapsed = true
      c_dom = Domain.new(rep_doms[0].from, prev_dom.to, prev_dom.did, -1, self.pid, prev_dom.cid, "#{repno} collapsed repeats")
      p.add_domain(c_dom)
      doms.push(c_dom)
    elsif (repno > 1)
      doms += rep_doms.push(prev_dom)
      rep_doms.push(prev_dom).each {|dom| p.add_domain(dom)}
    else
      doms.push(prev_dom)
      p.add_domain(prev_dom)
    end
#    return doms
    p.update_arrstr()
    return p
  end

  # collapse in place
  def collapse! (repunit=2)
    p = self.collapse(repunit)
    self.comment += "Repeats collapsed if found"
    @domains = p.domains
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

  protected

  def update_arrstr
    @arr_str = @domains.collect{|d| d.did}.join(@sep)
  end

  private

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

  def initialize (from, to, did, evalue, pid, clanid=nil, comment=String.new, acc=String.new, go=String.new, sequence=String.new)
    @from = from
    @to = to
    @did = did
    @acc = acc
    @go = go
    @cid = clanid
    @evalue = evalue
    @sequence = sequence
    #@protein = protein
    @pid = pid
    @type = type_check
    @comment = comment
    @length = @to - @from
  end
  attr_accessor :from, :to, :did, :acc, :go, :cid, :evalue, :sequence, :type, :pid, :comment, :length

  def fasta_seq
    return (sequence.empty?) ?
      false :
      ">#{self.pid}.#{self.did}\t#{self.from}-#{self.to}\n#{self.sequence}"
  end

  def to_s(clans = false)
    did = self.did
    if clans
      did = (self.cid.nil?) ? did : cid
    end
    self.from.to_s+"\t"+self.to.to_s+"\t"+did+"\t"+self.evalue.to_s
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
