This is FAPROTAX (Functional Annotation of Prokaryotic Taxa), a database that maps prokaryotic clades (e.g. genera, species or subspecies) to established metabolic or other ecologically relevant functions based on the current literature. FAPROTAX includes software for converting taxonomic microbial community profiles (e.g. in the form of an OTU table) into putative functional profiles, based on taxa identified in a sample.


What is included
- - - - - - - - -
1. This very informative README file.
2. The FAPROFAX database ('FAPROFAX.txt').
3. A python script ('collapse_table.py') for transforming taxonomic tables into functional tables based on the taxa shared with the FAPROTAX database.

Please read the license agreement at the end of this document prior to using FAPROTAX.


Caveats
- - - -
FAPROTAX defines functional groups in terms of taxa (e.g. species or genera) affiliated with each functional group. These affiliations are mostly based on peer-reviewed literature, such as announcements of cultured representatives. The authors of FAPROTAX made a great effort to ensure the accuracy of functional annotations. However, an implicit assumption of FAPROTAX is that if all cultured members of a clade can perform a particular function (at the time of the publication cited), then all members of the clade (cultured and non-cultured) can perform that function. As more organisms are being cultured, some of these generalization may turn out to be false.


Usage instructions
- - - - - - - - - -
To convert an OTU table (BIOM format) into a functional table, use the shell command:
   ./collapse_table.py -i otu_table.biom -o functional_table.biom -g FAPROTAX.txt
The above command assumes that observation IDs in the BIOM table are full taxon names (e.g. SILVA or Greengenes format). 

Alternatively, if observation IDs are OTU numbers and taxonomic annotation is instead included as observation metadata ('taxonomy'), use the following command:
   ./collapse_table.py -i otu_table.biom -o functional_table.biom -g FAPROTAX.txt --collapse_by_metadata 'taxonomy'

It is also possible to convert classical tables (e.g. in TSV format). An overview of all options can be obtained via the script's help menu, by calling:
   ./collapse_table.py -h

For more detailed usage instructions please visit the official FAPROTAX website at:
   http://www.zoology.ubc.ca/louca/FAPROTAX

If you need to modify the FAPROTAX database, for example to add new functional annotations, please consult the instructions and the license agreement within the database file.




LICENSE AGREEMENT
- - - - - - - - -
Copyright (c) 2020, Stilianos Louca
All rights reserved.

Use and redistributions of the FAPROTAX database, including any associated software, with or without modification, are permitted provided that the following conditions are met:

   * Redistributions must retain the above copyright notice, this list of conditions and the following disclaimer in the database itself, as well as in documentation and/or other materials provided with the distribution.

   * Neither the name of the original author (Stilianos Louca), nor the names of its contributors may be used to endorse or promote products derived from FAPROTAX without specific prior written permission.

   * If the FAPROTAX database or any associated software has been modified from its original version, this needs to be clearly indicated in any redistribution and in any publication using the modified version.

THE FAPROTAX DATABASE, INCLUDING ANY ASSOCIATED SOFTWARE, IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF FAPROTAX, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



