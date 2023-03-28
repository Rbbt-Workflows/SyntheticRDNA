Simulates reads from ribosomal RNA contigs

Uses (NEAT GenReads)[https://github.com/zstephens/neat-genreads] to simulate sequencing of ribosomal genes.
Currently it only simulates DNA. It can simulate the T2T sample or a fully synthetic sample in which the
ribosomal gene contigs are taken from a synthetic catalogue of morphs, each generated from one of the 24
unique morphs found in the T2T with a random set of mutations. The mutations used for the synthetic catalogue
for a randomly generated catalogue that is used by each synthetic morph to pick a random subset; this way
morphs can share mutations across them

# Tasks

## simulate_t2t

Simulates the t2t sample using a FASTA file with all the contigs

## align_simulated_t2t

Align the simulated t2t sample using the HTS workflow

## mutation_catalogue

Creates a catalogue of mutations that will be incorporated into the different synthetic morphs

The mutations are listed with padding around them, instead of by position, so that the same mutation may
be applied to slightly different morph templates.

## morph_catalogue

Generate the morph catalogue to be used in simulating the sample FASTA

Uses a collection of template morphs and generates a new morph catalogue by taking random morphs from the
reference catalogue and introducing randomly mutations from the mutation catalogue

## sample_fasta

Simulate a FASTA file 

The file consists of a number of contigs formed by taking morphs from the catalogue with replacement

## simulate_sample

Simulates sequencing reads from a sample FASTA

## simulate_sample_cohort

Simulates sequencing reads from a cohort of samples
