=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     https://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.
 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut


=head1 ImportRNAi

This module imports RNAi data.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportRNAi;

use strict;
use warnings;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source_info;

sub fetch_input {
  #create output folder structure and fetches input files
  my $self = shift;
  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $species      = $self->required_param('species');

  $self->debug($self->param('debug_mode'));
  $self->core_db_adaptor($self->get_species_adaptor('core'));
  $self->variation_db_adaptor($self->get_species_adaptor('variation'));

 # import specific constants

  %source_info = (source_description => 'An RNAi screen',
                  object_type => 'Gene',
		  somatic_status => 'somatic',
		  source_name => 'RNAi',
		  source_name_short => 'RNAi',
    		  source_version => '1'
		);

  #create workdir folder
  my $workdir = $pipeline_dir."/".$source_info{source_name_short}."/".$species;
  make_path($workdir);
  $self->workdir($workdir); 

  open (my $logFH, ">", $workdir."/".'log_import_out_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open (my $errFH, ">", $workdir."/".'log_import_err_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  open (my $pipelogFH, ">", $workdir."/".'log_import_debug_pipe_'.$source_info{source_name_short}.'_'.$species) || die ("Failed to open file: $!\n");
  $self->logFH($logFH);
  $self->errFH($errFH);
  $self->pipelogFH($pipelogFH);
	
  #find data file
  my $inputFile = 'RNAi_phenotypes.txt';
  die "expected $workdir/$inputFile" unless -e $workdir."/".$inputFile;

  $self->param('RNAi_file', $inputFile);

}

sub run {
  my $self = shift;
  my $input_file = $self->required_param('RNAi_file'); 

  # get phenotype data
  my $results = $self->parse_input_file($input_file);
  $self->print_logFH("Got ".(scalar @{$results->{'phenotypes'}})." phenotypes \n") if ($self->debug);

  # save phenotypes
  $self->save_phenotypes(\%source_info, $results);

  my %param_source = (source_name => $source_info{source_name},
                      type => $source_info{object_type});

  $self->param('output_ids', { source => \%param_source,
                               species => $self->required_param('species')
                             });
}


sub write_output {
  my $self = shift;

  $self->print_pipelogFH("Passing $source_info{source_name_short} import (".$self->required_param('species').") for checks (check_phenotypes)\n") if ($self->debug);
  close($self->logFH) if defined $self->logFH ;
  close($self->errFH) if defined $self->errFH ;
  close($self->pipelogFH) if defined $self->pipelogFH ;
  $self->dataflow_output_id($self->param('output_ids'), 1);
}


=head2 parse_input_file
  Arg [1]    : string $infile
               The input file name.
  Example    : $results = $obj->parse_input_file($infile)
  Description: Parse phenotypes from RNAi input file
  Returntype : hashref with results (key 'phenotypes')
  Exceptions : none

=cut


sub parse_input_file {
  my ($self, $infile) = @_;

  my @phenotypes;

  #Open the input file for reading
  open (IN,'<',$self->workdir."/".$infile) || die ("Could not open $infile for reading\n");

  #Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    my @data = split /\t/, $_;

    die "couldn't parse input file" unless scalar(@data) == 3;

  #For now - until we can infer gene from dsRNA alignment

    my $gene_id     	 = $data[0];
    my $phenotype   	 = $data[1];
    my $study	         = $data[2];

    push @phenotypes, {
    'id' 	  => $gene_id,
    'description' => $phenotype,
    'study'	  => $study,
    'external_id' => $gene_id
    };

  }
  close(IN);
  return {'phenotypes' => \@phenotypes };
}

1
  

