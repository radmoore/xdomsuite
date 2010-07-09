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
  def pfamscan(evalue=10, name=true, length=true)
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
      
      #next if (/^#{@comment}/.match(line) || /^$/.match(line))

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
class Proteome

  HEADERre = /^>(\S+)\s+(\d+)$/
  FASTAHEAD = /^>(\S+)\s?.?/
  include Enumerable


  # Returns a new xdom object containing proteins and domains
  #
  # * +filename+ : path to the xdomfile
  # * +evalue+   : annotation threshold
  # * +species+  : name of the proteome (OPTIONAL)
  def initialize(filename, evalue, species = "")
    @comment = '#'
    @filename = filename
    @file_ext = File.extname(@filename)
    @file_bn = File.basename(@filename, @file_ext)
    @species  = (species.empty?) ? get_species(filename) : species
    @evalue   = evalue
    @proteins = Hash.new
    @arrangements = Hash.new
    @domains  = Hash.new
    @total_domains    = 0
    @residue_coverage = 0.00
    @total_dom_residues = 0
    @total_prot_length = 0
    @clans = false
    @names = true
    if (@file_ext == '.hmmout')
      read_pfamscan(evalue)
    elsif (@file_ext == '.xdom')
      read_xdom(@filename)
    else
      raise "Unsupported file type #{@file_ext}"
    end
    @total_proteins   = @proteins.size
    @uniq_domains   = @domains.size
    @current_prot = 0
    @domains2go = nil
  end

  attr_reader :filename, :evalue, :species, :total_proteins, :total_domains, :uniq_domains, :clans, :names

  def clans=(var)
    return if var == @clans
    unless (var.kind_of?(TrueClass) or var.kind_of?(FalseClass))
      raise "Class variable 'clans' must of type boolean"
    end
    @clans = var
    @proteins.values.each {|p| p.clans = var}
    update_arrangements()
    update_domains()
  end

  def names=(var)
    return if var == @names
    unless (var.kind_of?(TrueClass) or var.kind_of?(FalseClass))
      raise "Class variable 'names' must of type boolean"
    end
    @names = var
    @proteins.values.each {|p| p.names = var}
    update_arrangements()
    update_domains()
  end

  # return a list of all domain ids that occur at least once in each proteome
  def intersection_of_domains(other_proteome)
    return @domains.keys & other_proteome.domains_hash.keys
  end

  # return a list of all domain ids that occur at least once in one of the proteomes
  def union_of_domains(other_proteome)
    return @domains.keys | other_proteome.domains_hash.keys
  end

  # returns a random protein object
  # or nil
  def grab
    return nil if @proteins.keys.empty?
    return @proteins[@proteins.keys[rand(@proteins.size)]]
  end

  # returns all proteins that have more than 1 domain
  def get_multidomain_proteins
    return @proteins.values.select{|p| p.is_multidomain?}
  end

  # Returns string representation of self
  def to_s
    lines = String.new
    @proteins.values.each {|p| 
      p.clans = @clans
      p.names = @names
      lines += p.to_s
    }
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

  def domains_hash
    return @domains
  end


  # Returns an array of domain identifiers representing the unique set of domains present in self 
  def get_all_uniq_doms
    @domains.keys
  end

  # wrapper around same method in
  # protein
  def pfam_A
    doms = Array.new
    @proteins.values.each{|p| doms += p.pfam_A }
    return doms
  end

  # wrapper around same method in
  # protein
  def pfam_B
    doms = Array.new
    @proteins.values.each{|p| doms += p.pfam_B }
    return doms
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

  # TODO: this does not work - seems to modify self (no deep copy?)!!!
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
    while(line = f.gets)
      line.chomp!
      if (m = FASTAHEAD.match(line))
        pid = m[1]
        self.add_protein(Protein.new(pid)) unless @proteins.member?(pid)
      else
        seq = line
        @proteins[pid].sequence += seq
      end
    end
    f.close
  end

  # Provides an interface to attach an annot file from BLAST2GO which allows to map
  # GO terms to domains via the proteins
  #
  # Will add the possiblity to collect all GO terms of a domain from all proteins it
  # occurs in.
  def attach_blast2go_annot(file)
    f = File.open(file, "r")
    while(line = f.gets)
      line.chomp!
      pid, goid = line.split("\t")[0,2]
      next unless @proteins.include?(pid)
      @proteins[pid].goterms[goid] = 1 unless @proteins[pid].goterms.include?(goid)
    end
    update_domains2go()
  end

  def update_domains2go
    @domains2go = Hash.new([])
    @domains.each do |did, pids|
      pids.each{|pid| @domains2go[did] += @proteins[pid].goterms.keys}
    end
  end

  def domain2go(did)
    return @domains2go[did]
  end

  # Provides an interface to attach a GFF file to the Domain annotation in _self_
  #
  # Requires that gene annotation is reported by maker, i.e. source (2nd column) == "maker".
  # Pfam domains are then added as attributes to the mRNA entry that codes for the Protein
  # containing the individual domain.
  def attach_gff3(file)
    f = File.open(file, "r")
    while(line = f.gets)
      line.chomp!
      next if line[0,1] == "#"
      next unless line.split("\t").count == 9
      seqid, source, type, start, stop, score, strand, phase, attributes = line.split("\t")
      next unless source == "maker"
      cl = ChromosomalLocation.new(seqid, source, type, start, stop, score, strand, phase, attributes)
      next unless cl.type == "mRNA"
      attrs = cl.attributes_hash
      next unless attrs['ID'] or attrs['Parent']
      next unless @proteins[attrs['ID']] or @proteins[attrs['Parent']]
      proteins2check = Array.new
      proteins2check << @proteins[attrs['ID']] if attrs['ID'] and @proteins.include?(attrs['ID'])
      proteins2check << @proteins[attrs['Parent']] if attrs['Parent'] and @proteins.include?(attrs['Parent'])
      proteins2check.each do |p|
        p.domains.each{|d| cl.attributes_string += "Pfam=#{d.acc};"}
        p.chromosomal_locations[cl.type] = cl
      end
    end
  end

  # outputs domain annotation of the proteome in gff3 format
  def to_gff3
    @proteins.each do |pid, protein|
      protein.chromosomal_locations.each{|type, cl| puts cl.to_s}
    end
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
      next if uniq.has_key?(p.arrstr)
      uniq[p.arrstr] = p
    end
    return uniq.values
  end

  # Returns the procentage of amino acid residues that fall into an annotated region (domain). 
  # If _num_ is _true_ , a float is returned. Otherwise a formatted string is returned.
  #
  # TODO:
  # * check number return
  def res_coverage(num = false)
    cov = (@total_dom_residues.to_f/@total_prot_length.to_f)
    return num ? cov : sprintf('%.2f', cov)
  end

  # Returns the procentage of proteins that contain at least one domain. 
  # If _num_ is _true_ , a float is returned. Otherwise a formatted string is returned.
  #
  # TODO:
  # * check number return
  def prot_coverage(num = false)
    protwdoms = (@proteins.values.select{|p| p.has_domains?}).count
    cov = ( protwdoms.to_f / @total_proteins.to_f )
    return num ? cov : sprintf('%2f', cov)
  end

  def resolve_overlaps(mode)
    @proteins.values.each {|p| p.resolve_overlaps(mode)}
		update_arrangements()
    return nil
  end

  def simple_overlap_resolution
    total_domains = 0
    @proteins.values.each {|p| 
      p.simple_overlap_resolution
      total_domains += p.domains.length
    }
    update_domains()
		update_arrangements()
    @total_domains = total_domains
    return nil
  end

  def resolve_overlaps_with_sets
    @proteins.values.each {|p| p.resolve_overlaps_with_sets}
    update_arrangements()
    return nil
  end


  # TODO
  def annotate_with_context(xdom)

  end

  def add_protein(p)
    @proteins[p.pid] = p unless @proteins.member?(p.pid)
    @total_proteins += 1
  end


private

  # TODO: perhaps this can update domains as well
	def update_arrangements
		new_arr = Hash.new
		@proteins.each do |pid, p|
			new_arr[p.arrstr] = Array.new unless (new_arr.has_key?(p.arrstr))
			new_arr[p.arrstr].push(pid)
		end
		@arrangements = new_arr
	end

  # TODO: when clan or domains is changed in an instance of
  # proteome, recreate domain hash - ensure that this cannot be
  # done in update_arrangements()
  # this will be painfully slow due to use of includes?
  # * consider using uniq instead of includes?
  def update_domains
    new_doms = Hash.new
    @proteins.each do |pid, p|
      p.domains.each do |d|
        new_doms[d.did] = Array.new unless (new_doms.has_key?(d.did))
        new_doms[d.did] << pid unless new_doms[d.did].include?(pid)
      end
    end
    @domains = new_doms
    @uniq_domains = @domains.size
    return true
  end

  def get_species(filename)
    file = File.basename(filename, ".xdom")
    file.split('.')[0]
  end

  # TODO:
  # * !!!!!sanitize
  # * check efficiency !!!
  def read_pfamscan(evalue=10, envelope=false)
    hmmout = File.open(@filename, "r")
    cpid = pid = nil
    entries = Array.new
    while(line = hmmout.gets)
      next if (/^#{@comment}/.match(line) || /^$/.match(line))
      line.chomp!
      cpid = line.split[1]
      eva_ht = line.split[13]
      next if (eva_ht.to_f > evalue.to_f)
      if (pid.nil?)
        entries << line
        pid = cpid
        next
      end
      if (cpid != pid)
        if entries.empty?
          entries = Array.new
          entries << line
          pid = cpid
          next
        end
        p = nil
        entries.each do |domline|
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
           pre_as) = domline.split  # predicted_active_site_residues
          hmm_ac = hmm_ac.split('.')[0] if (/.+\.\d+/.match(hmm_ac))
          cla_id = nil if (cla_id == 'NA' || cla_id == 'No_clan')
          if p.nil?
            p = Protein.new(seq_id, seq_le.to_i, "", @species, "", @filename)
            p.clans = @clans
            p.names = @names
          end
          from = (envelope) ? env_st.to_i : aln_st.to_i
          to = (envelope) ? env_en.to_i : aln_en.to_i         
          d = Domain.new(from, to, hmm_na, eva_ht.to_f, p.pid, cla_id, "", hmm_ac)
          p.add_domain(d)
          @total_dom_residues += (to - from)
          did = (@names) ? hmm_na : hmm_ac
          @domains[did] = Array.new unless (@domains.member?(did))
          @total_domains += 1 if (@domains[did] << p.pid)
        end
        @arrangements[p.arrstr] = Array.new unless (@arrangements.has_key?(p.arrstr))
        @arrangements[p.arrstr] << p.pid
        @total_prot_length += p.length
        @proteins[p.pid] = p
        entries = Array.new
        entries << line
        pid = cpid
        next
      end
      entries << line
      pid = cpid
    end
    p = nil
    entries.each do |domline|
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
       pre_as) = domline.split  # predicted_active_site_residues
      hmm_ac = hmm_ac.split('.')[0] if (/.+\.\d+/.match(hmm_ac))
      cla_id = nil if (cla_id == 'NA' || cla_id == 'No_clan')
      if p.nil?
        p = Protein.new(seq_id, seq_le.to_i, "", @species, "", @filename)
        p.clans = @clans
        p.names = @names
      end

      p = Protein.new(seq_id, seq_le.to_i, "", @species, "", @filename) if p.nil?
      from = (envelope) ? env_st.to_i : aln_st.to_i
      to = (envelope) ? env_en.to_i : aln_en.to_i         
      d = Domain.new(from, to, hmm_na, eva_ht.to_f, p.pid, cla_id, "", hmm_ac)
      p.add_domain(d)
      @total_dom_residues += (to - from)
      did = (@names) ? hmm_na : hmm_ac
      @domains[did] = Array.new unless (@domains.member?(did))
      @total_domains += 1 if (@domains[did] << p.pid)
    end

    @arrangements[p.arrstr] = Array.new unless (@arrangements.has_key?(p.arrstr))
    @arrangements[p.arrstr] << p.pid
    @total_prot_length += p.length
    @proteins[p.pid] = p

    hmmout.close

  end

  def read_xdom (xdomfile)
    f = File.open(xdomfile, "r")
    protein = 0
      while(line = f.gets)
        line.chomp!
        next if (/^$/.match(line) || /^#/.match(line))
        if (m = HEADERre.match(line)) then
          unless (protein == 0)
            @arrangements[protein.arrstr] = Array.new unless (@arrangements.has_key?(protein.arrstr))
            @arrangements[protein.arrstr].push(protein.pid)
          end
          protein = Protein.new(m[1], m[2].to_i, "", @species, "", @filename)
          @total_prot_length += m[2].to_i
          @proteins[m[1]] = protein
          next
        end
        (from, to, did, evalue) = line.split
        clan = nil
        domain = Domain.new(from.to_i, to.to_i, did, evalue.to_f, protein.pid, clan)
        @total_dom_residues += (to.to_i-from.to_i)
        @domains[did] = Array.new unless (@domains.member?(did))
        @total_domains += 1 if (@domains[did].push(protein.pid))
        protein.add_domain(domain)
      end
      @arrangements[protein.arrstr] = Array.new unless (@arrangements.has_key?(protein.arrstr))
      @arrangements[protein.arrstr].push(protein.pid)
    f.close
  end


end

# CLASS: CHROMOSOMALLOCATION
class ChromosomalLocation
  def initialize(seqid, source, type, start, stop, score, strand, phase, attributes)
    @seqid = seqid
    @source = seqid
    @type = type
    @start = start.to_i
    @stop = stop.to_i
    @score = score.to_f
    @strand = strand
    @phase = phase
    @attributes_string = attributes
    update_attributes_hash(@attributes_string)
  end

  attr_reader :seqid, :source, :type, :start, :stop, :score, :strand, :phase, :attributes_string, :attributes_hash

  def attributes_string=(var)
    @attributes_string = var
    update_attributes_hash(var)
  end

  def attributes_hash=(var)
    @attributes_hash = var
    update_attributes_string(var)
  end

  def to_s
    return [@seqid, @source, @type, @start.to_s, @stop.to_s, @score.to_s, @strand, @phase, @attributes_string].join("\t")
  end

  private

  def update_attributes_hash(attr_string)
    @attributes_hash = Hash.new
    attr_string.split(";").each do |a|
      key, value = a.split("=")
      @attributes_hash[key] = value
    end
  end

  def update_attributes_string(attr_hash)
    @attributes_string = ""
    keys = attr_hash.keys
    @attributes_string << "ID=#{attr_hash['ID']};" if attr_hash.include?('ID')
    @attributes_string << "Parent=#{attr_hash['Parent']};" if attr_hash.include?('Parent')
    keys.each do |k|
      next if k == 'ID' or k == 'Parent'
      @attributes_string << "#{k}=#{attr_hash[k]};" 
    end
  end

end

# CLASS: PROTEIN
# TODO
# * attr_reader instead of accessor
class Protein

  FASTAHEAD = /^>(\S+)\s+.+/

  def initialize (pid, length = 0, sequence=String.new, species=String.new, comment=String.new, srcfile=String.new)
    @pid      = pid
    @length   = length
  	@sequence = sequence
    @species  = species
    @comment  = comment
    @srcfile  = srcfile
    @domains  = Array.new
    @arrstr  = String.new
    @current_dom = 0
    @position = Hash.new
		@deleted = nil
		@collapsed = false
    @sep = ';'
    @clans = true
    @names = true
    @pfamb = false
    @chromosomal_locations = Hash.new
    @goterms = Hash.new
  end
  attr_accessor :length, :species, :comment, :sep, :pfamb, :chromosomal_locations, :goterms
	attr_reader :deleted, :pid, :arrstr, :clans, :names, :sequence
  
  def is_multidomain?
    return @domains.size > 1
  end

  def add_domain (domain)
    @domains.push(domain)
    self.update_arrstr()
  end

  def clans=(var)
    return if (var == @clans)
    unless (var.kind_of?(TrueClass) or var.kind_of?(FalseClass))
      raise "Class variable 'clans' must of type boolean"
    end
    @clans = var
    self.update_arrstr() if @clans != var
  end

  def names=(var)
    return if (var == @names)
    unless (var.kind_of?(TrueClass) or var.kind_of?(FalseClass))
      raise "Class variable 'names' must of type boolean"
    end
    @names = var
    self.update_arrstr()
  end

  def domains
    return @domains
  end

  # get a random domain
  def grab
    return nil if @domains.empty?
    return @domains[rand(@domains.size)]
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
    return @domains.select{|d| d.type == 'A'}
		#p = self.type_filter('A')
		#return p.domains
	end

  # Wrapper for type_filter. Returns an array of Domain objects of type PfamB, or an empty array if no PfamB domains are found
	def pfam_B
    return @domains.select{|d| d.type == 'B'}
		#p = self.type_filter('B')
		#return p.domains
	end

  # Returns a copy of _self_, containing only Domains of _type_. If no Domain objects of _type_ are present in _self_,
  # then self.domains will return an empty array
	def type_filter(type)
		p = Protein.new(@pid, @length, @sequence, @species, @comment)
		@domains.each {|d| p.add_domain(d) if (d.type == type)}
		p.update_arrstr()
		return p
	end

  # In place variant of self.type_filter
	def type_filter!(type)
		doms_filtered = Array.new
		@domains.each {|d| doms_filtered << d if (d.type == type) }
		@domains = doms_filtered
		self.update_arrstr()
		return self
	end

  def member?(did)
    @domains.each {|d| return true if (d.did == did)}
    return false
  end

	def collapsed?
		return @collapsed
	end

  def has_overlaps?
    cdom = pdom = nil
    pos = 0
    until (pos == self.total_domains)
      cdom = self.domains[pos]
      if (pdom.nil?)
        pdom = cdom
        pos += 1
        next
      end
      return true if (pdom.to >= cdom.from)
      pos += 1
    end
    return false
  end

  # Create all non-overlapping domain sets and find 
  # the set that maximized coverage
  def resolve_overlaps_with_sets
    return self if self.total_domains < 2 # overlaps impossible
  #  puts self.pid
    v = true if (self.pid == 'BGIBMGA005287-PA')
    non_overlapping = Array.new
    pos = 0
    until (pos == self.total_domains) 
      non_overlapping << [self.domains[pos]]
      non_overlapping << find_domainsets(self.domains[pos], pos)
      pos += 1
    end

   # if (self.pid == 'BGIBMGA005287-PA')
   #   c = 0
   #   non_overlapping.each do |pos|
   #     puts "SET: #{c}"
   #     non_overlapping[c].each do |d|
   #       puts "\t"+d.to_s
   #     end
   #     puts
   #     c += 1
   #   end
   # end
   # exit


    # determine best set
    maxpos = coverage = pos = 0
    until (pos == non_overlapping.size)
      c = 0
      domarray = non_overlapping[pos]
 #     puts domarray.inspect if v
      domarray.each {|d| c += d.length}
      if (c > coverage)
        coverage = c
        maxpos = pos
      end
      pos += 1
    end
    if (self.pid == 'BGIBMGA005287-PA')
      c = 0
      non_overlapping.each do |pos|
        puts "SET: #{c}"
        non_overlapping[c].each do |d|
          puts "\t"+d.to_s
        end
        puts
        c += 1
      end
    end
    @domains = non_overlapping[maxpos]
    self.update_arrstr
    return self
  end

  # find non-overlapping domain sets
  def find_domainsets(cdom, pos)
    v = true if (self.pid == 'BGIBMGA005287-PA')
    puts "POS: #{pos}, cdom: #{cdom}" if v
    cpos = 0
    noover = Array.new
    noover << cdom
    until (cpos == self.total_domains)
      if cpos == pos  
        cpos += 1
        next
      end
      dom = self.domains[cpos]
      #if (cdom.overlaps?(dom))
      if (cpos < pos)
        if (cdom.to >= dom.from)
          puts "#{cdom.did} [#{cdom.to}] overlaps with #{dom.did} [#{dom.from}]" if v
        else
          noover << dom
        end
      else
        if (cdom.from >= dom.to)
          puts "#{cdom.did} [#{cdom.to}] overlaps with #{dom.did} [#{dom.from}]" if v
        else
          noover << dom
        end
      end
      cpos += 1
    end
    return (noover.empty?) ? nil : noover
  end

#  def overlaps?(domain)
#    raise "overlap? requires formal paramter of type <DOMAIN>" unless (domain.instance_of?(Domain))
#    return true if (self.to >= domain.from)
#    return false
#  end  


  def simple_overlap_resolution
    pos = 0
    pdom = cdom = nil
    deleted = Array.new
    until (pos == self.total_domains)
      cdom = self.domains[pos]
      if (pdom.nil?)
        pdom = cdom
        pos += 1
        next
      end
      if (pdom.to > cdom.from)
        if (pdom.type == cdom.type)
          if (cdom.evalue <= pdom.evalue)
            self.domains.delete(pdom)
            deleted << pdom
          else
            self.domains.delete(cdom)
            deleted << cdom
          end
        else 
          if (pdom.type == 'B')
            self.domains.delete(pdom)
            deleted << pdom
          elsif (cdom.type == 'B')
            self.domains.delete(cdom)
            deleted << cdom
          else
            if (cdom.evalue <= pdom.evalue)
              self.domains.delete(pdom)
              deleted << pdom
            else
              self.domains.delete(cdom)
              deleted << cdom
            end
          end
        end
        cdom = pdom = nil
        pos = 0
        next
      end
      pos += 1
      pdom = cdom
    end
    self.update_arrstr
    return self
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
		self.update_arrstr
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
      did = (@names) ? d.did : d.acc
      did = (@clans && (not d.cid.nil?)) ? d.cid : did
      doms += "#{d.from.to_s}\t#{d.to.to_s}\t#{did}\t#{d.evalue.to_s}"
      doms += (d.comment.empty?) ? "\n" : "\t;#{d.comment}\n"
    end
    head+doms
  end

  # TODO: check if required (why not attr_reader?)
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

  def collapse(repunit=2)
    #v = true if self.pid == 'gi|66523667|ref|XP_625264.1|'
    return self.dup if self.total_domains < repunit
    repno = 1
    prev_dom = nil
    rep_doms = Array.new
   # doms = Array.new
    p = Protein.new(@pid, @sequence, @species, @comment)
    self.domains.each do |d|

      unless (prev_dom.nil?)
        if (prev_dom.did == d.did)
          repno += 1
          rep_doms.push(prev_dom)
        else
          if (repno >= repunit)
						@collapsed = true

            c_dom = Domain.new(rep_doms[0].from, prev_dom.to, prev_dom.did, -1, self.pid, prev_dom.cid, "#{repno} collapsed repeats", prev_dom.acc)
            p.add_domain(c_dom)
           # doms.push(c_dom)
            rep_doms.clear
          elsif (repno > 1)
           # doms += rep_doms.push(prev_dom)
            rep_doms.push(prev_dom).each {|dom| p.add_domain(dom)}
            rep_doms.clear
          else
            p.add_domain(prev_dom)
           # doms.push(prev_dom)
          end
          repno = 1
        end
      end
      prev_dom = d
    end
    if (repno >= repunit)
			@collapsed = true
      c_dom = Domain.new(rep_doms[0].from, prev_dom.to, prev_dom.did, -1, self.pid, prev_dom.cid, "#{repno} collapsed repeats", prev_dom.acc)
      p.add_domain(c_dom)
     # doms.push(c_dom)
    elsif (repno > 1)
      #doms += rep_doms.push(prev_dom)
      rep_doms.push(prev_dom).each {|dom| p.add_domain(dom)}
    else
      #doms.push(prev_dom)
      p.add_domain(prev_dom)
    end
    p.update_arrstr()
    return p
  end

  # collapse in place
  def collapse! (repunit=2)
    p = self.collapse(repunit)
    self.comment += "Repeats collapsed if found (REPUNIT #{repunit})"
    @domains = p.domains
    self.update_arrstr()
    return
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
    arr1 = self.arrstr.split(';').sort
    arr2 = protein.arrstr.split(';').sort
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

  def jaccard_dist(protein, collapse=true)
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

  def get_interdoms(min=20)
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

  def sequence=(seq)
    @sequence = seq
    @length = seq.length
    @domains.each{|d| d.sequence = seq[(d.from-1)..(d.to-1)]}
  end

  protected

  def update_arrstr
    return if @domains.nil?
    @arrstr = @domains.collect{|d|
      setid = (@names) ? d.did : d.acc
      setid = (@clans && (not d.cid.nil?)) ? d.cid : setid
      setid
    }.join(@sep)
  end

  private

  #### PUT SET FINDER HERE####

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

  ATYPE_ACC = /PF.+/
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

  # Compares end position of_self_ to start position of _domain_.
  # Returns true if positions overlap, false otherwise
  def overlaps?(domain)
    raise "overlap? requires formal paramter of type <DOMAIN>" unless (domain.instance_of?(Domain))
    return true if (self.to >= domain.from)
    return false
  end

  private

  def type_check
    return 'A' if (ATYPE_ACC.match(self.acc))
    return 'B' if (BTYPE_ID.match(self.did))
    return 'B' if (BTYPE_ACC.match(self.acc))
    return 'C' if (CTYPE.match(self.did))
    return 'D' if (DTYPE.match(self.did))
    return 'A'
  end


end
