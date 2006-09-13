#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::IndividualAdaptor

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  my $ia = $db->get_IndividualAdaptor();

  # Get an individual by its internal identifier
  my $ind = $ia->fetch_by_dbID(52);

  # Get all individuals with a particular name
  foreach my $ind (@{$ia->fetch_all_by_name('PKH053(M)')}) {
    print "Individual ", $ind->name(), "\n";
  }

  # get all individuals from a population
  my $pop = $pop_adaptor->fetch_by_name('PACIFIC');
  foreach my $ind (@{$ia->fetch_all_by_Population($pop)}) {
    print $ind->name(), "\n";
  }

  # get all children of an individual
  foreach my $child (@{$ia->fetch_all_by_parent($ind)}) {
    print $child->name(), " is a child of ", $ind->name(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for Individual objects.
Individuals may be retrieved from the ensembl variation database by
several means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor;

use Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::Individual;
our @ISA = ('Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor');


=head2 fetch_individual_by_synonym

    Arg [1]              : $individual_synonym
    Example              : my $ind = $ind_adaptor->fetch_individual_by_synonym($individual_synonym,$source);
    Description          : Retrieves individual for the synonym given in the source. If no source is provided, retrieves all the synonyms
    Returntype           : list of Bio::EnsEMBL::Variation::Individual
    Exceptions           : none
    Caller               : general

=cut

sub fetch_individual_by_synonym{
    my $self = shift;
    my $synonym_name = shift;
    my $source = shift;
    my $individuals;
    my $ind;
    #return all sample_id from the database
    my $samples = $self->SUPER::fetch_sample_by_synonym($synonym_name, $source);
    foreach my $sample_id (@{$samples}){
	#get the ones that are individuals
	$ind = $self->fetch_by_dbID($sample_id);
	push @{$individuals}, $ind if (defined $ind);
    }
    return $individuals;
}

=head2 fetch_all_by_name

  Arg [1]    : string $name the name of the individuals to retrieve
  Example    : my @inds = @{$ind_adaptor->fetch_all_by_name('CEPH1332.05')};
  Description: Retrieves all individuals with a specified name.  Individual
               names may be non-unique which is why this method returns a
               reference to a list.
  Returntype : reference to a list of Individual ids
  Exceptions : throw if no argument passed
  Caller     : general

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = shift;

  defined($name) || throw("name argument expected");

  my $sth = $self->prepare
    (q{SELECT i.sample_id, s.name, s.description,
              i.gender, i.father_individual_sample_id, i.mother_individual_sample_id, it.name, it.description
       FROM   individual i, sample s, individual_type it
       WHERE  s.name = ?
       AND    it.individual_type_id = i.individual_type_id
       AND    s.sample_id = i.sample_id});

  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->execute();

  my $result =  $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}




=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : my $pop = $pop_adaptor->fetch_by_name('PACIFIC');
               foreach my $ind (@{$ia->fetch_all_by_Population($pop)}) {
                 print $ind->name(), "\n";
               }
  Description: Retrieves all individuals from a specified population 
  Returntype : reference to list of Bio::EnsEMBL::Variation::Individual objects
  Exceptions : throw if incorrect argument is passed
               warning if provided Population does not have an dbID
  Caller     : general

=cut

sub fetch_all_by_Population {
  my $self = shift;
  my $pop = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw("Bio::EnsEMBL::Variation::Population arg expected");
  }

  if(!$pop->dbID()) {
    warning("Population does not have dbID, cannot retrieve Individuals");
    return [];
  }

  my $sth = $self->prepare
    (q{SELECT i.sample_id, s.name, s.description,
              i.gender, i.father_individual_sample_id, i.mother_individual_sample_id, it.name, it.description
       FROM   individual i, individual_population ip, sample s, individual_type it
       WHERE  i.sample_id = ip.individual_sample_id
       AND    i.sample_id = s.sample_id
       AND    i.individual_type_id = it.individual_type_id
       AND    ip.population_sample_id = ?});

  $sth->bind_param(1,$pop->dbID,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;
}




=head2 fetch_all_by_parent_Individual

  Arg [1]    : Bio::EnsEMBL::Variation::Individual
  Example    : my @children = @{$ia->fetch_all_by_parent_Individual($ind)};
  Description: Retrieves all individuals which are children of a provided
               parent individual.  This function operates under the assumptions
               that Male individuals can only be fathers, Female individuals
               can only be mothers and Unknown individuals can only be one
               or the other - not both.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Individuals
  Exceptions : throw if incorrect argument passed
               warning if provided individual has no dbID 
  Caller     : general, Individual::get_all_child_Individuals

=cut

sub fetch_all_by_parent_Individual {
  my $self = shift;
  my $parent  = shift;

  if(!ref($parent) || !$parent->isa('Bio::EnsEMBL::Variation::Individual')) {
    throw("Bio::EnsEMBL::Variation::Individual argument expected");
  }

  if(!defined($parent->dbID())) {
    warning("Cannot fetch child Individuals for parent without dbID");
    return [];
  }

  my $gender = $parent->gender() || '';

  my $father_sql =
    q{SELECT i.sample_id, s.name, s.description,
             i.gender, i.father_individual_sample_id, i.mother_individual_sample_id, it.name, it.description
      FROM   individual i, sample s, individual_type it
      WHERE  i.father_individual_sample_id = ?
      AND    i.individual_type_id = it.individual_type_id
      AND    s.sample_id = i.sample_id};
  my $mother_sql =
    q{SELECT i.sample_id, s.name, s.description,
              i.gender, i.father_individual_sample_id, i.mother_individual_sample_id, it.name, it.description
       FROM   individual i, sample s, individual_type it
       WHERE  i.mother_individual_sample_id = ?
       AND    i.individual_type_id = it.individual_type_id
       AND    i.sample_id = s.sample_id};

  if($gender eq 'Male') {
    my $sth = $self->prepare($father_sql);
    $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
    $sth->execute();
    my $result = $self->_objs_from_sth($sth);
    $sth->finish();
    return $result;
  }
  elsif($gender eq 'Female') {
    my $sth = $self->prepare($mother_sql);
    $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
    $sth->execute();
    my $result = $self->_objs_from_sth($sth);
    $sth->finish();
    return $result;
  }

  # unknown gender

  my $sth = $self->prepare($mother_sql);
  $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
  $sth->execute();
  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  # if this parent was a mother, finish now and return results
  return if(@$result);

  # otherwise assume was a father (or nothing)
  $sth = $self->prepare($father_sql);
  $sth->bind_param(1,$parent->dbID,SQL_INTEGER);
  $sth->execute();
  $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return $result;
}

=head2 fetch_all_strains

    Args       : none
    Example    : my $strains = $ind_adaptor->fetch_all_strains();
    Description: Retrieves Individuals that should be considered as strain (fully inbred) in the specie.
    Returntype : list of Bio::EnsEMBL::Variation::Individual
    Exceptions : none
    Caller     : Bio:EnsEMBL:Variation::Individual

=cut

sub fetch_all_strains{
    my $self = shift;
    
    return $self->generic_fetch("it.name = 'fully_inbred'");

}

=head2 get_display_strains

    Args       : none
    Example    : my $strains = $ind_adaptor->get_display_strains();
    Description: Retrieves strain_names that are going to be displayed in the web (reference + default + others)
    Returntype : list of strings
    Exceptions : none
    Caller     : web

=cut

sub get_display_strains{
    my $self = shift;
    my @strain_names;
    my $name;
    #first, get the reference strain
    $name = $self->get_reference_strain_name();
    push @strain_names, $name;
    #then, get the default ones
    my $default_strains = $self->get_default_strains();
    push @strain_names, @{$default_strains};
    #and finally, get the others
    my $sth = $self->prepare(qq{SELECT meta_value from meta where meta_key = ?
				});
    $sth->bind_param(1,'individual.display_strain',SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$name);
    while ($sth->fetch()){
	push @strain_names, $name;
    }
    $sth->finish;
    return \@strain_names;

}


=head2 get_default_strains

    Args       : none
    Example    : my $strains = $ind_adaptor->get_default_strains();
    Description: Retrieves strain_names that are defined as default in the database(mainly, for web purposes)
    Returntype : list of strings
    Exceptions : none
    Caller     : web

=cut

sub get_default_strains{
    my $self = shift;
    my @strain_names;
    my $name;
    my $sth = $self->prepare(qq{SELECT meta_value from meta where meta_key = ?
				});
    $sth->bind_param(1,'individual.default_strain',SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$name);
    while ($sth->fetch()){
	push @strain_names, $name;
    }
    $sth->finish;
    return \@strain_names;

}


=head2 get_reference_strain_name

    Args       : none
    Example    : my $reference_strain = $ind_adaptor->get_reference_strain_name();
    Description: Retrieves the reference strain_name that is defined as default in the database(mainly, for web purposes)
    Returntype : string
    Exceptions : none
    Caller     : web

=cut

sub get_reference_strain_name{
    my $self = shift;

    my $name;
    my $sth = $self->prepare(qq{SELECT meta_value from meta where meta_key = ?
				});
    $sth->bind_param(1,'individual.reference_strain',SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$name);
    $sth->fetch();
    $sth->finish;

    return $name;

}

#
# private method, constructs Individuals from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($dbID, $name, $desc, $gender, $father_id, $mother_id,$it_name,$it_desc);

  $sth->bind_columns(\$dbID, \$name, \$desc, \$gender,
                     \$father_id, \$mother_id, \$it_name, \$it_desc);

  my %seen;
  my %wanted_fathers;
  my %wanted_mothers;

  my @inds;


  while($sth->fetch()) {
    # get objects for mother and father if they were already constructed
    # otherwise may have to be lazy-loaded later

    my $father;
    if(defined($father_id)) {
      $father = $seen{$father_id};
      if(!$father) {
        $wanted_fathers{$father_id} ||= [];
        push @{$wanted_fathers{$father_id}}, $father_id;
      }
    }
    my $mother;
    if(defined($mother_id)) {
      $mother = $seen{$mother_id};
      if(!$mother) {
        push @{$wanted_mothers{$mother_id}}, $mother_id;
      }
    }


    my $ind = $seen{$dbID} ||= Bio::EnsEMBL::Variation::Individual->new
      (-dbID        => $dbID,
       -adaptor     => $self,
       -description => $desc,
       -gender      => $gender,
       -name        => $name,
       -father_individual => $father,
       -mother_individual => $mother,
       -father_individual_id => $father_id,
       -mother_individual_id => $mother_id,
       -type_individual => $it_name,
       -type_description => $it_desc);

    $seen{$dbID} = $ind;

    push @inds, $ind;
  }

  # load any of the 'wanted' parent individuals that we did not have at the
  # of creation, but which we have now

  foreach my $wanted_id (keys %wanted_fathers) {
    if($seen{$wanted_id}) {
      # add father to every child that wanted it
      foreach my $ind_id (@{$wanted_fathers{$wanted_id}}) {
        $seen{$ind_id}->father_Individual($seen{$wanted_id});
      }
    }
  }
  foreach my $wanted_id (keys %wanted_mothers) {
    if($seen{$wanted_id}) {
      # add mother to every child that wanted it
      foreach my $ind_id (@{$wanted_mothers{$wanted_id}}) {
        $seen{$ind_id}->mother_Individual($seen{$wanted_id});
      }
    }
  }

  return \@inds;

}

sub _tables{return (['individual','i'],
		    ['sample','s'],
		    ['individual_type','it'])}

sub _columns{
    return qw(s.sample_id s.name s.description i.gender i.father_individual_sample_id i.mother_individual_sample_id it.name it.description);
}

sub _default_where_clause{
    return 's.sample_id = i.sample_id AND i.individual_type_id = it.individual_type_id';
}

1;
