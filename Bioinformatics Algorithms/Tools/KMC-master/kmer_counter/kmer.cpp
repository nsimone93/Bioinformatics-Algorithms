#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.0.0
  Date   : 2017-01-28
*/

#include "defs.h"
#include "kmer.h"

uint32 CKmer<1>::QUALITY_SIZE      = 0;
uint32 CKmerQuake<1>::QUALITY_SIZE = 4;

uint32 CKmer<1>::KMER_SIZE = 1;
uint32 CKmerQuake<1>::KMER_SIZE = 1;
// ***** EOF
