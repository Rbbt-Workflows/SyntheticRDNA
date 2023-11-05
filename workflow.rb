require 'rbbt-util'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/SyntheticRDNA'

Workflow.require_workflow "NEATGenReads"
module SyntheticRDNA
  extend Workflow

  helper :load_morphs do |fasta_file|
    morphs = {}
    name_line = nil
    TSV.traverse fasta_file, :type => :array do |line|
      line.strip!
      if line.chomp.empty?
        next
      elsif line.start_with? ">"
        name_line = line
        morphs[name_line] = ""
      else
        morphs[name_line] += line.strip
      end
    end

    morphs
  end

  helper :preprocess_mutation do |morph,pad,regions,size=1|
    pass = false
    while ! pass
      pos = rand(morph.length).floor
      pass = true
      regions.each do |region|
        if region.include? pos
          pass = false
        elsif (pos..pos+size-1).include?(region.first)
          pass = false
        end
      end if regions
    end

    start = pos - pad
    eend = pos + pad
    start = 0 if start < 0

    reference = morph[pos] 
    original = morph[start..eend] 

    [pos, start, eend, reference, original]
  end

  input :reference_morphs_for_mutations, :file, "FASTA with reference morphs", Rbbt.data["T2T_rDNA.consensus.single.fa.gz"].find
  input :number_of_snvs, :integer, "Number of SNV to introduce", 100
  input :number_of_ins, :integer, "Number of insertions", 20
  input :number_of_del, :integer, "Number of deletions", 20
  input :pad, :integer, "Surronding area", 10
  input :conserved_regions, :text, "BED file with regions to conserve"
  task :mutation_catalogue => :array do |reference_morphs,number_of_snvs,number_of_ins,number_of_del,pad,conserved_regions|
    morphs = load_morphs reference_morphs

    regions = {}
    conserved_regions.split("\n").each do |line| 
      morph_code, s,e = line.split("\t")
      region = (s.to_i..e.to_i) 
      regions[morph_code] ||= []
      regions[morph_code] << region
    end if conserved_regions

    snvs = number_of_snvs.times.collect do 
      morph_code = morphs.keys.sample
      morph = morphs[morph_code]
      pos, start, eend, reference, original = preprocess_mutation morph, pad, regions[morph_code[1..-1]]

      alt = (%w(A C T G) - [reference.upcase]).shuffle.first
      alt = alt.downcase if reference == reference.downcase
      mutated = begin
                  tmp = morph.dup
                  tmp[pos] = alt
                  tmp[start..eend] 
                end
      [original, mutated] * "=>" + "-" + [morph_code, pos + 1, [reference,alt]*">",pos - start] * ":"
    end

    ins = number_of_ins.times.collect do 
      morph_code = morphs.keys.sample
      morph = morphs[morph_code]
      pos, start, eend, reference, original = preprocess_mutation morph, pad, regions[morph_code]

      size = rand(6).to_i + 2
      #pos = [pos, morph.length - 1].min  # Added by aisha Ensure pos is within the valid range

      alt = size.times.collect{ %w(A C T G).sample } * ""
      ref_allele = morph[pos] # added by aisha
      mutated = begin
                  tmp = morph.dup
                  tmp[pos] += alt
                  tmp[start..eend+size] 
                end
      #[original, mutated] * "=>" + "-" + [morph_code, pos + 1, "+" + alt,pos - start] * ":"
      #[morph[pos], mutated] * "=>" + "-" + [morph_code, pos, "+" + alt, pos - start] * ":" # added by aisha
      [original, mutated] * "=>" + "-" + [morph_code, pos + 1, ref_allele + ">" + ref_allele + alt,pos - start] * ":"  # added by aisha
    end

    dels = number_of_del.times.collect do 
      morph_code = morphs.keys.sample
      morph = morphs[morph_code]

      size = rand(6).to_i + 1
      pos, start, eend, reference, original = preprocess_mutation morph, pad, regions[morph_code], size
      #pos = [pos, morph.length - size].min  # Added by aisha Ensure pos is within the valid range
      mutated = begin
                  tmp = morph.dup
		  deleted_sequence = tmp[(pos-1)..(pos+size-1)]
                  ref_allele = deleted_sequence[0]
		  #deleted_sequence = tmp[(pos..pos+size-1)] # added by aisha
		  #deleted_sequence = tmp[(pos-1..pos+size-1)] # added by aisha
                  tmp[(pos..pos+size-1)] = ""
                  tmp[start..(eend-size)] 
                end
      #[original, mutated] * "=>" + "-" + [morph_code, pos + 1, "-" * size,pos - start] * ":"
      [original, mutated] * "=>" + "-" + [morph_code, pos , deleted_sequence +">"+ ref_allele , pos - start] * ":" # added by aisha
    end

    snvs + ins + dels
  end

  dep :mutation_catalogue
  input :reference_morphs_for_templates, :file, "FASTA with reference morphs", Rbbt.data["T2T_rDNA45S.24_uniq_morphs.fa.gz"].find
  input :catalogue_size, :integer, "Number of synthetic morphs to create", 144
  input :mutations_per_morph, :integer, "Number of mutations to introduce in each morph", 10
  extension "fa"
  task :morph_catalogue => :text do |reference_morphs,catalogue_size,mutations_per_morph|
    original_morphs = load_morphs reference_morphs
    mutations = step(:mutation_catalogue).load.collect{|e| e.split(/=>|->/) }

    pad = step(:mutation_catalogue).inputs[:pad]

    original_morph_keys = original_morphs.keys
    catalogue_size.times.collect do |morph_number|
      morph_source = original_morph_keys.sample
      original_sequence = original_morphs[morph_source]

      name = ">synth_#{morph_number}.#{morph_source[1..-1]}"

      selected_mutations = mutations
        .select{|ref,mut| original_sequence.include? ref }
        .sample(mutations_per_morph)
        .collect{|ref,mut,info| [ref, mut, info, original_sequence.index(ref)] }
        .sort_by{|ref,mut,info,pos| pos }

      mutation_info = selected_mutations.collect{|ref,mut,info| 
        orig_morph, orig_pos, change = info.split(":")
        orig_pos = orig_pos.to_i

        fixed_pad = [pad, orig_pos].min
        location = original_sequence.index(ref) + 1 + fixed_pad
        fixed = [morph_source[1..-1], location, change] * ":"

        offset = selected_mutations
          .collect{|ref,mut,info| info }
          .reject{|i| i.split(":")[2].include? ">" }
          .select{|i| i.split(":")[1].to_i < orig_pos }
          .inject(0) do |acc,i|  
            c = i.split(":")[2]
            shift = c.include?("+") ? c.length - 1 : - c.length
            acc += shift
          end

        fixed_and_offset = [name[1..-1], location + offset, change] * ":"


        [info, location, fixed, fixed_and_offset] * "\t" 
      }

      mutated_sequence = original_sequence.dup
      acc_offset = 0
      selected_mutations.each do |ref,mut,info,pos|
        offset = info.split(":").last.to_i

        mutated_sequence[pos+offset+acc_offset..pos+acc_offset+ref.length-1] = mut[offset..-1]

        acc_offset += mut.size - ref.size
      end

      file('mutations')[name[1..-1]].write(mutation_info * "\n" + "\n")

      [name, mutated_sequence] * "\n"
    end * "\n"
  end

  dep :morph_catalogue, :jobname => "SharedCatalogue"
  input :sample_contigs, :integer, "Number of sample contigs", 200
  extension "fa.gz"
  task :sample_fasta_balanced => :text do |sample_contigs|
    catalogue = load_morphs step(:morph_catalogue).path
    catalogue_keys = catalogue.keys
    selected_base_morphs = []
    txt = sample_contigs.times.collect do |sample_contig|
      base = catalogue_keys.sample
      selected_base_morphs << base
      sequence = catalogue[base]
      name = ">sample_#{sample_contig}.#{base[1..-1]}"
      [name, sequence] * "\n"
    end * "\n"
    tmpfile = file('tmp.fa')
    Open.write(file('selected_base_morphs.list'), selected_base_morphs * "\n")
    Open.write(tmpfile, txt + "\n")
    CMD.cmd("bgzip #{tmpfile}")
    Open.mv tmpfile + '.gz', self.tmp_path
    nil
  end

  dep :morph_catalogue, :jobname => "SharedCatalogue"
  input :sample_contigs_min, :integer, "Min number of contigs in sample", 200
  input :sample_contigs_max, :integer, "Max number of contigs in sample", 600
  input :number_of_base_morphs_min, :integer, "Min number of different base morphs in sample", 5
  input :number_of_base_morphs_max, :integer, "Max number of different base morphs in sample", 20
  input :number_of_copies_per_base_morph_min, :integer, "Min number of copies per base morph", 1
  input :number_of_copies_per_base_morph_max, :integer, "Max number of copies per base morph", 100
  extension "fa.gz"
  task :sample_fasta => :text do |sample_contigs_min,sample_contigs_max, number_of_base_morphs_min, number_of_base_morphs_max, number_of_copies_per_base_morph_min, number_of_copies_per_base_morph_max|
    catalogue = load_morphs step(:morph_catalogue).path
    catalogue_keys = catalogue.keys

    sample_contigs = (sample_contigs_min..sample_contigs_max).to_a.sample
    number_of_base_morphs = (number_of_base_morphs_min..number_of_base_morphs_max).to_a.sample

    Log.medium "sample_contigs: #{sample_contigs}"

    base_morphs = []
    while base_morphs.length < number_of_base_morphs
      base_morphs << catalogue_keys.sample
      base_morphs.uniq!
    end

    selected_base_morphs = []
    sample_contig = 1
    txts = []
    while selected_base_morphs.length < sample_contigs
      base_morphs.each_with_index do |base,base_num|
        copies = (number_of_copies_per_base_morph_min..number_of_copies_per_base_morph_max).to_a.sample
        Log.medium "base_morph_copies #{[base, copies] * ": "}"
        copies.times do |copy|
          selected_base_morphs << base
          sequence = catalogue[base]
          name = ">sample_#{sample_contig}.#{base[1..-1]}"
          txts << [name, sequence] * "\n"
          sample_contig += 1
        end
        break if selected_base_morphs.length > sample_contigs_max
        next if base_num + 1 < number_of_base_morphs_min 
        break if selected_base_morphs.length > sample_contigs
      end
    end

    Log.medium "sample_contigs real: #{selected_base_morphs.length}"

    tmpfile = file('tmp.fa')
    Open.write(file('selected_base_morphs.json'), Misc.counts(selected_base_morphs).to_json)
    Open.write(file('selected_base_morphs.list'), selected_base_morphs * "\n" + "\n")
    Open.write(tmpfile, txts * "\n" + "\n")
    CMD.cmd("bgzip #{tmpfile}")
    Open.mv tmpfile + '.gz', self.tmp_path
    nil
  end


  dep :sample_fasta
  dep_task :simulate_sample, NEATGenReads, :NEAT_simulate_DNA, :reference => :sample_fasta

  input :number_of_samples, :integer, "How many samples to generate", 100
  dep :simulate_sample do |jobname,options|
    options[:number_of_samples].to_i.times.collect do |i|
      {:inputs => options, :jobname => "Sample#{i}"}
    end
  end
  task :simulate_sample_cohort => :array do
    dependencies.collect{|dep| dep.load }.flatten
  end

  input :number_of_samples, :integer, "How many samples to generate", 100
  dep :sample_fasta do |jobname,options|
    options[:number_of_samples].to_i.times.collect do |i|
      {:inputs => options, :jobname => "Sample#{i}"}
    end
  end
  task :sample_cohort_fasta => :array do
    dependencies.collect{|dep| dep.path }.flatten
  end

end

#require 'SyntheticRDNA/tasks/basic.rb'

#require 'rbbt/knowledge_base/SyntheticRDNA'
#require 'rbbt/entity/SyntheticRDNA'

