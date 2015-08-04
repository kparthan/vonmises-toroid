#include "pdbClass.h"
#include <unistd.h>

/**********************************************************************
*                       A Brookhaven PDB format.                      *
*---------------------------------------------------------------------*
*                     FORMAT FOR ATOM RECORD                          * 
* --------------------------------------------------------------------*
* COL   ARR_IND  DATA TYPE     FIELD        DEFINITION                *
* --------------------------------------------------------------------*
*  1-6  | 0-5   | Record name |  "ATOM "  |                           *
*  7-11 | 6-10  | Integer     |  serial   | Atom serial number.       *
* 13-16 | 12-15 | Atom        |  name     | Atom name.                *
* 17    | 16    | Character   |  altLoc   | Alt. location indicator.  *
* 18-20 | 17-19 | Residue name|  resName  | Residue name.             *
* 22    | 21    | Character   |  chainID  | Chain identifier.         *
* 23-26 | 22-25 | Integer     |  resSeq   | Residue sequence number.  *
* 27    | 26    | AChar       |  iCode    | Code for res insertion.   *
* 31-38 | 30-37 | Real(8.3)   |  x        | Coords for X in Angstroms.*
* 39-46 | 38-45 | Real(8.3)   |  y        | Coords for Y in Angstroms.*
* 47-54 | 46-53 | Real(8.3)   |  z        | Coords for Z in Angstroms.*
* 55-60 | 54-59 | Real(6.2)   | occupancy | Occupancy.                *
* 61-66 | 60-65 | Real(6.2)   | tempFactor| Temperature factor.       *
* 77-78 | 76-77 | LString(2)  | element   | Elem symb,right-justif.   *
* 79-80 | 78-79 | LString(2)  | charge    | Charge on the atom.       *
*                                                                     *
***********************************************************************/

// empty global ATOM_t object
ATOM_t tmpatm ;


/*#### ATOM_t CLASS ###*/

ATOM_t::ATOM_t() {
	EMPTY = true ;
}

const double* ATOM_t::getXYZ() const {
	return &coord[0] ;
}

void ATOM_t::process( const char line[] ) {
	EMPTY = false ;

	//store line verbatim
	verbatimLine.assign( line ) ;
	string &verb = verbatimLine ;

	//store record type -- ATOM or HETATM?
	recordType = verb.substr( 0, 6 ) ;
	
	//store atom serial number
	string tmp ;
	tmp = verb.substr( 6, 5 ) ;
	serial = (uint32_t)( atoi(tmp.c_str()) ) ;

	//store atom name
	atomName = verb.substr(12, 4 ) ;

	//store altLoc
	altLoc = line[16] ;

	//store resName
	resiName = verb.substr( 17, 3 ) ;

	//store chainID
	chainID = line[21] ;

	//store resi sequence number
	tmp = verb.substr( 22, 4 ) ;
	resiSeqNum = (int32_t)( atoi(tmp.c_str()) ) ;

	//store icode
	iCode = line[26] ;

	//store x,y,z coords
	//x-coord
	tmp = verb.substr( 30, 8 ) ;
	sscanf( tmp.c_str() , "%lf", &x) ;
	coord[0] = x ;
	//y-coord
	tmp = verb.substr( 38, 8 ) ;
	sscanf( tmp.c_str() , "%lf", &y) ;
	coord[1] = y ;
	//z-coord
	tmp = verb.substr( 46, 8 ) ;
	sscanf( tmp.c_str() , "%lf", &z) ;
	coord[2] = z ;

	//store occupancy
	tmp = verb.substr( 54, 6 ) ;
	sscanf( tmp.c_str() , "%lf", &occupancy) ;

	//store temp factor
	tmp = verb.substr( 60, 6 ) ;
	sscanf( tmp.c_str() , "%lf", &tempFactor) ;
}

bool ATOM_t::isEmpty() {
	return EMPTY ;
}

string ATOM_t::getResIdentifier() const {
	stringstream tmp ;
	tmp  << resiName << resiSeqNum  ;
	return tmp.str() ;
}

int32_t ATOM_t::getResSeqNum() const {
	return resiSeqNum ; 
}

double* ATOM_t::getCoords() {
	return coord ;
}

void ATOM_t::print() {
	cout << endl ;
	cout  << "\t\t\t\t" << "isempty? " << isEmpty() << endl ;
	cout << verbatimLine << endl ;
	cout <<  "\t\t\t\t" << "record Type " << recordType << endl ;
	cout  << "\t\t\t\t" << "serial number " << serial << endl ;
	cout  << "\t\t\t\t" << "atomname " << atomName << endl ;
	cout  << "\t\t\t\t" << "altLoc " << altLoc << endl ;
	cout  << "\t\t\t\t" << "resiName " << resiName << endl ;
	cout  << "\t\t\t\t" << "chain ID " << chainID << endl ;
	cout  << "\t\t\t\t" << "resiSeqNum " << resiSeqNum << endl ;
	cout  << "\t\t\t\t" << "icode " << iCode << endl ;
	cout  << "\t\t\t\t" << "coords " << std::fixed
		<< coord[0] << " " 
		<< coord[1] << " " 
		<< coord[2] << endl ;
	cout  << "\t\t\t\t" << "Occupancy " << occupancy << endl ;
	cout  << "\t\t\t\t" << "tempFactor " << tempFactor << endl ;
}

void ATOM_t::printPDB(ofstream &out) {
  out << verbatimLine << endl;
}

void ATOM_t::printPDB(ofstream &out, Superpose3DClass &supobj) {
  vector<double> v;
  v.push_back(coord[0]); v.push_back(coord[1]); v.push_back(coord[2]);
  supobj.transformVector(v);

  out << verbatimLine.substr(0,30);
  out << fixed << setw(8) << setprecision(3) << v[0]; 
  out << fixed << setw(8) << setprecision(3) << v[1]; 
  out << fixed << setw(8) << setprecision(3) << v[2]; 
  out << verbatimLine.substr(54,26) << endl;
}



/*#### RESIDUE_t CLASS ###*/
RESIDUE_t::RESIDUE_t() {
	nAtoms = 0 ;
	ssType = '-' ; // default value indicating that ss assignment 
		       //is incomplete.
	ssDefFlag = false ;

	dihedralsPhiPsiOmega[0] = dihedralsPhiPsiOmega[1] 
		= dihedralsPhiPsiOmega[2] = INFINITY ;

	NHCO_HBondResIndx[0] = NHCO_HBondResIndx[1] = ~(0u) ;
}

const uint32_t RESIDUE_t::getnAtoms() const {
	return nAtoms ;
}

/* Member function to get an atom record.
 * On success, an atom object identified by the pattern is returned.
 * On failure, that is when the pattern is invalid, a NULL object is returned.
 * 
 * Note that the pattern is searched for a substring match among the  
 * keys of the atomName_indx map. (Keys are of type string.) 
 * So if one were to search for a C_ALPHA atom, the pattern string
 * can be "CA" or " CA" or " CA ". Note that the latter two are more
 * specific that the former. Use patterns with caution. */
const ATOM_t& RESIDUE_t::getAtomRecord( const char pattern[] ) {
	map<string,int>::iterator it ;

	bool foundflag =false ;
	for ( it=atomName_indx.begin() ; it != atomName_indx.end(); it++ ) {
		size_t stat = (*it).first.find( pattern ) ;
		if (stat != string::npos ) {
			foundflag = true ;
			return atom[ (*it).second ] ;
		}	
	}

	// At this stage, the pattern string is not a part of the 
	// any atoms record. 
	// Return a reference to an empty global ATOM_t object.
	return tmpatm ;
}

/* Member function to get coords of a specified atom.
 * On success, a pointer to double[3] corresponding to xyz coordinates is returned.
 * On failure, that is when the pattern is invalid, a null pointer is returned.
 * 
 * Note that the pattern is searched for a substring match among the  
 * keys of the atomName_indx map. (Keys are of type string.) 
 * So if one were to search for a C_ALPHA atom, the pattern string
 * can be "CA" or " CA" or " CA ". Note that the latter two are more
 * specific that the former. Use patterns with caution. */
double* RESIDUE_t::getCoords( const char pattern[] ) {
	map<string,int>::iterator it ;

	bool foundflag =false ;
	for ( it=atomName_indx.begin() ; it != atomName_indx.end(); it++ ) {
		size_t stat = (*it).first.find(pattern) ;
		if (stat != string::npos ) {
			foundflag = true ;
			return atom[(*it).second].getCoords() ;
		}	
	}
	return NULL ;
}


/* Member function to return CA coords.
 * NOTE: On failure, that is when no CA record exists, 
 * the returns NULL.
 * */
double* RESIDUE_t::getCACoords() {
	map<string,int>::iterator it ;

	bool foundflag =false ;
	for ( it=atomName_indx.begin() ; it != atomName_indx.end(); it++ ) {
		size_t stat = (*it).first.find("CA") ;
		if (stat != string::npos ) {
			foundflag = true ;
			return atom[(*it).second].getCoords() ;
		}	
	}
	return NULL ;
}


char RESIDUE_t::getResCode() const {
	return resCode ;
}


string RESIDUE_t::getResIdentifier() const {
	stringstream tmp ;
	if(nAtoms > 0 ) {
		tmp << '(' << atom[0].getResIdentifier() << ')' ;
	}
	return tmp.str() ;
}

int32_t RESIDUE_t::getResSeqNum() const {
	if(nAtoms > 0 ) {
		return atom[0].getResSeqNum() ; 
	}
	else return ~0 ; 
}


void RESIDUE_t::aaAssign( const char aa[], char code ) {
	strcpy( resType, aa ) ;
	resCode = code ;
}

bool RESIDUE_t::isAtomIndexed( const string atomName ) {
	static map<string,int>::iterator itr;
	itr = atomName_indx.find( atomName ) ;
	if( itr == atomName_indx.end() ) return false ;
	
	else return true ;
}

void RESIDUE_t::process( const char line[] ) {
	static char recordstr[7] ;
	static char atmstr[6] ;

	strncpy( recordstr, line, 6 ) ; recordstr[6] = '\0' ;

	atmstr[0] = line[12] ;
	atmstr[1] = line[13] ;
	atmstr[2] = line[14] ;
	atmstr[3] = line[15] ;
	atmstr[4] = line[16] ;
	atmstr[5] = '\0'     ;

	if( isAtomIndexed( atmstr ) == false ) {
		//insert into map, atomstr-->vectorindx
		atomName_indx.insert(pair<string,int>(atmstr,nAtoms));
		nAtoms++ ;
		//append a new chain object
		ATOM_t newatom ;
		atom.push_back( newatom ) ;
		
		int currAtomIndx = nAtoms-1 ;
		atom[currAtomIndx].process( line ) ;
	}
	//else if the line contains an atom type that is already indexed
	// in this atom object, then ignore it. Such ambiguous repeats
	//  are in non-conformity  with the Brookhaven format.
}

void RESIDUE_t::setPhiPsiOmega(
  const double phi,
  const double psi,
  const double omega
) {
	dihedralsPhiPsiOmega[0] = phi ;
	dihedralsPhiPsiOmega[1] = psi ;
	dihedralsPhiPsiOmega[2] = omega ;
}

void RESIDUE_t::setNHCO_HBondResIndx( const uint32_t resIndx, const uint32_t arrIndx) {
	assert( arrIndx < 2 ) ;
	NHCO_HBondResIndx[arrIndx] = resIndx ;
}

uint32_t* RESIDUE_t::getNHCO_HBondResIndx() {
	return NHCO_HBondResIndx ;
}


double* RESIDUE_t::getPhiPsiOmega() {
	return dihedralsPhiPsiOmega ;
}

void RESIDUE_t::print() {
	cout << endl ;
	cout  << "\t\t\t" << "atom_name indx: \n" ;
	map<string,int>::iterator it ;
	cout << "\t\t\t" << "RES TYPE: " << resType ;
	cout << "(" << resCode << ")" << endl ;

	//ss related
	if( ssDefFlag ) {
		cout  << "\t\t\t" << "sec. str. type " << ssType << endl ;
	}

	for ( it=atomName_indx.begin() ; it != atomName_indx.end(); it++ )
		cout  << "\t\t\t" << (*it).first << " => " 
			<< (*it).second << endl;
	cout << endl ;

	cout  << "\t\t\t" << "nAtoms = " << nAtoms << endl ;
	for( unsigned int i = 0 ; i < nAtoms ; i++ ) {
		cout  << "\t\t\t" << "ATOM " << i << endl ;
		atom[i].print() ;
	}
}

void RESIDUE_t::printPDB(ofstream &out) {
	for( unsigned int i = 0 ; i < nAtoms ; i++ ) {
		atom[i].printPDB(out) ;
	}
}

void RESIDUE_t::printPDB(ofstream &out, Superpose3DClass &supobj) {
	for( unsigned int i = 0 ; i < nAtoms ; i++ ) {
		atom[i].printPDB(out,supobj) ;
	}
}



// secondary structure related functions
char RESIDUE_t::getSSType() const {
	assert(ssDefFlag) ;
	return ssType ;
}

void RESIDUE_t::setSSType( const char c ) {
	ssDefFlag = true ;
	ssType = c;
}

/*#### CHAIN_t CLASS ###*/
CHAIN_t::CHAIN_t() {
	nResidues = 0 ;
	chainID = (char )0 ;
	aacode["ALA"] = 'A'; 
	aacode["ARG"] = 'R'; 
	aacode["ASN"] = 'N';
	aacode["ASP"] = 'D'; 
	aacode["ASX"] = 'B'; 
	aacode["CYS"] = 'C';
	aacode["SEC"] = 'U';
	aacode["GLU"] = 'E'; 
	aacode["GLN"] = 'Q'; 
	aacode["GLX"] = 'Z'; 
	aacode["PCA"] = '!'; //my own code for rare pyroglutamic acid
	aacode["GLY"] = 'G';
	aacode["HIS"] = 'H'; 
	aacode["ILE"] = 'I';
	aacode["LEU"] = 'L'; 
	aacode["XLE"] = 'J';
	aacode["LYS"] = 'K'; 
	aacode["MET"] = 'M';
	aacode["MSE"] = '$'; //my own code for selenomethionine
	aacode["PHE"] = 'F'; 
	aacode["PRO"] = 'P';
	aacode["SER"] = 'S'; 
	aacode["THR"] = 'T'; 
	aacode["TRP"] = 'W';
	aacode["TYR"] = 'Y'; 
	aacode["VAL"] = 'V'; 

	aacode["PYL"] = 'O'; 
	aacode["XAA"] = 'X';		
	aacode["UNK"] = 'X'; // alternative to UNK		

	//some strange amino acid 3 letter codes. 
	//Treat them as unknown
	aacode["INI"] = 'X'; // alternative to UNK		
	aacode["HYP"] = 'X'; // alternative to UNK		
	
	//Some proteins in the PDB are Protein-DNA complexes
	// and are crystallized with protein bound to the DNA
	// ATOM records of the DNA part have 
	// " DA", " DC", " DG", " DT" in the amino acid column. 
	// Treat them with special symbol '#'
	aacode[" DA"] = '#'; // my own code for DNA nucleotide, 'A'		
	aacode[" DT"] = '#'; // my own code for DNA nucleotide, 'T'		
	aacode[" DG"] = '#'; // my own code for DNA nucleotide, 'G'		
	aacode[" DC"] = '#'; // my own code for DNA nucleotide, 'C'		
	//some strange base in DNA
	aacode[" DI"] = '#'; // my own code for DNA		

	//Some PDBs have DNA/RNA bases represented as 
	// "  A", "  T", "  U", "  G", "  C"
	aacode["  A"] = '#'; // my own code for DNA/RNA nucleotide, 'A'		
	aacode["  T"] = '#'; // my own code for DNA/RNA nucleotide, 'T'		
	aacode["  U"] = '#'; // my own code for DNA/RNA nucleotide, 'U'		
	aacode["  G"] = '#'; // my own code for DNA/RNA nucleotide, 'G'		
	aacode["  C"] = '#'; // my own code for DNA/RNA nucleotide, 'C'		
	//some strange bases in DNA
	aacode["  I"] = '#'; // my own code for DNA		
	aacode["  N"] = '#'; // my own code for DNA		

	// SSE related.
	nSSEs = 0 ;
	sseFlag = axesAndCMFlag = tableauFlag = false ;
	setPymolColors() ;

}

void CHAIN_t::setPymolColors() {
	pymolColor.push_back("cyan");
	pymolColor.push_back("firebrick");
	pymolColor.push_back("white");
	pymolColor.push_back("marine");
	pymolColor.push_back("tv_blue");
	pymolColor.push_back("greencyan");
	pymolColor.push_back("warmpink");
	pymolColor.push_back("smudge");
	pymolColor.push_back("br6");
	pymolColor.push_back("ruby");
	pymolColor.push_back("teal");
	pymolColor.push_back("hotpink");
	pymolColor.push_back("salmon");
	pymolColor.push_back("palecyan");
	pymolColor.push_back("dash");
	pymolColor.push_back("br0");
	pymolColor.push_back("lightpink");
	pymolColor.push_back("density");
	pymolColor.push_back("palegreen");
	pymolColor.push_back("lightblue");
	pymolColor.push_back("lime");
	pymolColor.push_back("nitrogen");
	pymolColor.push_back("gray");
	pymolColor.push_back("tv_yellow");
	pymolColor.push_back("skyblue");
	pymolColor.push_back("deeppurple");
	pymolColor.push_back("red");
	pymolColor.push_back("violet");
	pymolColor.push_back("limegreen");
	pymolColor.push_back("olive");
	pymolColor.push_back("wheat");
	pymolColor.push_back("raspberry");
	pymolColor.push_back("sand");
	pymolColor.push_back("paleyellow");
	pymolColor.push_back("deepolive");
	pymolColor.push_back("tv_green");
	pymolColor.push_back("deepsalmon");
	pymolColor.push_back("orange");
	pymolColor.push_back("dirtyviolet");
	pymolColor.push_back("limon");
	pymolColor.push_back("deepblue");
	pymolColor.push_back("br1");
	pymolColor.push_back("brown");
	pymolColor.push_back("magenta");
	pymolColor.push_back("violetpurple");
	pymolColor.push_back("pink");
	pymolColor.push_back("hydrogen");
	pymolColor.push_back("tv_orange");
	pymolColor.push_back("deepteal");
	pymolColor.push_back("yelloworange");
	pymolColor.push_back("oxygen");
	pymolColor.push_back("br8");
	pymolColor.push_back("sulfur");
	pymolColor.push_back("darksalmon");
	pymolColor.push_back("chartreuse");
	pymolColor.push_back("grey");
	pymolColor.push_back("deepsalmon");
	pymolColor.push_back("yellow");
	pymolColor.push_back("blue");
	//pymolColor.push_back("green"); // we use this color for loops 
	pymolColor.push_back("br4");
	pymolColor.push_back("br2");
	pymolColor.push_back("bluewhite");
	pymolColor.push_back("tv_red");
	pymolColor.push_back("br7");
	pymolColor.push_back("aquamarine");
	pymolColor.push_back("purple");
	pymolColor.push_back("lightteal");
	pymolColor.push_back("lightmagenta");
	pymolColor.push_back("black");
	pymolColor.push_back("brightorange");
	pymolColor.push_back("slate");
	pymolColor.push_back("chocolate");
	pymolColor.push_back("purpleblue");
	pymolColor.push_back("br3");
	pymolColor.push_back("carbon");
	pymolColor.push_back("forest");
	pymolColor.push_back("splitpea");
	pymolColor.push_back("br9");
	pymolColor.push_back("lightorange");
	pymolColor.push_back("br5");
	nPymolColors = pymolColor.size() ;
}

const uint32_t CHAIN_t::getnResidues() const{
	return nResidues ;
}

const uint32_t CHAIN_t::getnAtoms( const uint32_t resIndx ) const{
	assert( resIndx < nResidues ) ; 
	return residue[resIndx].getnAtoms() ; 
}

char CHAIN_t::getResCode( const uint32_t resIndx ) const {
	assert( resIndx < nResidues ) ; 
	return residue[resIndx].getResCode() ; 

}

string CHAIN_t::getResIdentifier( const uint32_t resIndx ) const {
	assert( resIndx < nResidues ) ; 
	return residue[resIndx].getResIdentifier() ; 
}

int32_t CHAIN_t::getResSeqNum( 
  const uint32_t resIndx 
) const 
{
	assert( resIndx < nResidues ) ; 
	return residue[resIndx].getResSeqNum() ; 
}

const ATOM_t& CHAIN_t::getAtomRecord( const uint32_t resIndx, const char pattern[] ) {
	assert( resIndx < nResidues ) ; 
	return residue[resIndx].getAtomRecord( pattern ) ; 
}

double* CHAIN_t::getCACoords( const uint32_t resIndx ) {
	assert( resIndx < nResidues ) ; 
	return residue[resIndx].getCACoords() ; 
}

char CHAIN_t::getChainID() {
	return chainID ;
}

bool CHAIN_t::isResidueIndexed( const char resnumstr[] ) {
	static map<string,int>::iterator itr;
	itr = resNum_indx.find( resnumstr ) ;
	if( itr == resNum_indx.end() ) {
		return false ;
	}
	else return true ;
}

string CHAIN_t::getAminoAcidSequenceOfChain() {
  string aaseq;
  for (uint32_t i = 0 ; i < nResidues; i++) {
    aaseq.push_back(residue[i].getResCode());
  }
  return aaseq;
}

vector<string> CHAIN_t::getResIdentifiersOfChain() {
  vector<string> resids;
  for (uint32_t i = 0 ; i < nResidues; i++) {
    resids.push_back(residue[i].getResIdentifier());
  }
  return resids;
}

vector<vector<double> > CHAIN_t::getCACoordsOfChain() {
  vector<vector<double> > CA;
  for (uint32_t i = 0 ; i < nResidues; i++) {
    string patt = "CA";
    ATOM_t tatm = residue[i].getAtomRecord(patt.c_str());
    if( tatm.isEmpty() == true ) continue ; 
    const double *tca = tatm.getXYZ() ;

    vector<double> tmp ;
    tmp.push_back( tca[0] ) ;
    tmp.push_back( tca[1] ) ;
    tmp.push_back( tca[2] ) ;

    CA.push_back( tmp ) ;
  }
  return CA;
}

void CHAIN_t::getCACoordsOfChain(
  vector<vector<double> > &CA,
  string &aaseq,
  vector<string> &resIDStrings,
  vector<uint32_t> &resObjectIDs
) {
  string patt = "CA";
  for (uint32_t i = 0 ; i < nResidues; i++) {
    ATOM_t tatm = residue[i].getAtomRecord(patt.c_str());
    if( tatm.isEmpty() == true ) continue ; 
    const double *tca = tatm.getXYZ() ;

    vector<double> tmp ;
    tmp.push_back( tca[0] ) ;
    tmp.push_back( tca[1] ) ;
    tmp.push_back( tca[2] ) ;

    CA.push_back( tmp ) ;

    aaseq.push_back(residue[i].getResCode());
    resIDStrings.push_back((char)chainID+residue[i].getResIdentifier());
    resObjectIDs.push_back(i);
  }
}

void CHAIN_t::computePhiPsiOmega() {
	if( nResidues ==0 ) return ;
	for( uint32_t i = 0 ; i < nResidues ; i++ ) {
		double phi, psi, omega ;
		double *N_im1, *CA_im1, *CO_im1; // xyz of i-1th res backbone atoms 
		double *N_i, *CA_i, *CO_i ;      // xyz of ith res backbone atoms
		double *N_ip1, *CA_ip1, *CO_ip1 ;// xyz of i+1th res backbone atoms
		//init to NULL
		N_im1 = CA_im1 = CO_im1 = NULL ;
		N_i = CA_i = CO_i = NULL ;
		N_ip1 = CA_ip1 = CO_ip1 = NULL ;


		 N_i = residue[i].getCoords(" N") ;
		CA_i = residue[i].getCoords(" CA") ;
		CO_i = residue[i].getCoords(" C ") ;


		if( i > 0 ) {
			 N_im1 = residue[i-1].getCoords(" N ") ;
			CA_im1 = residue[i-1].getCoords(" CA") ;
			CO_im1 = residue[i-1].getCoords(" C ") ;
		}

		if( i < nResidues-1 ) {
			 N_ip1 = residue[i+1].getCoords(" N ") ;
			CA_ip1 = residue[i+1].getCoords(" CA") ;
			CO_ip1 = residue[i+1].getCoords(" C ") ;
		}


		/*
		cout << "i= " << i << " " ;
		if( N_im1 != NULL && CA_im1 != NULL && CO_im1 != NULL ) { 
			cout << "N_im1 "  <<  N_im1[0] << " " <<  N_im1[1] << " " <<  N_im1[2] << endl ; 
			cout << "CA_im1 " << CA_im1[0] << " " << CA_im1[1] << " " << CA_im1[2] << endl ; 
			cout << "CO_im1 " << CO_im1[0] << " " << CO_im1[1] << " " << CO_im1[2] << endl ; 
		}
		else {
			cout << "im1 is NULL\n" ;
		}

		if( N_i != NULL && CA_i != NULL && CO_i != NULL ) { 
			cout << "N_i " << N_i[0] << " " << N_i[1] << " " << N_i[2] << endl ; 
			cout << "CA_i " << CA_i[0] << " " << CA_i[1] << " " << CA_i[2] << endl ; 
			cout << "CO_i " << CO_i[0] << " " << CO_i[1] << " " << CO_i[2] << endl ; 
		}
		else {
			cout << "i is NULL" ;
		}

		if( N_ip1 != NULL && CA_ip1 != NULL && CO_ip1 != NULL ) { 
			cout << "N_ip1 " << N_ip1[0] << " " << N_ip1[1] << " " << N_ip1[2] << endl ; 
			cout << "CA_ip1 " << CA_ip1[0] << " " << CA_ip1[1] << " " << CA_ip1[2] << endl ; 
			cout << "CO_ip1 " << CO_ip1[0] << " " << CO_ip1[1] << " " << CO_ip1[2] << endl ; 
		}
		else {
			cout << "ip1 is NULL\n" ;
		}
		//sleep(1) ;
		*/

		phi = psi = omega = INFINITY ;

		if(i == 0 ) {
			if( N_i != NULL && CA_i != NULL && CO_i != NULL 
				&& N_ip1 != NULL ) {
				computeDihedralAngle(N_i,CA_i,CO_i,N_ip1,psi);
			}
			
			if(CA_i != NULL && CO_i != NULL && N_ip1 != NULL
			 	&& CA_ip1 != NULL ) {
				computeDihedralAngle(CA_i,CO_i,N_ip1,CA_ip1,omega);
			}
		}
		else if( i == nResidues-1 ) {
			if( CO_im1 != NULL && N_i != NULL && CA_i != NULL 
				&& CO_i != NULL ) {
				computeDihedralAngle(CO_im1,N_i,CA_i, CO_i,phi);
			}
		}
		else {
			if( CO_im1 != NULL && N_i != NULL && CA_i != NULL 
				&& CO_i != NULL ) {
				computeDihedralAngle(CO_im1,N_i,CA_i, CO_i,phi);
			}
			if( N_i != NULL && CA_i != NULL && CO_i != NULL 
				&& N_ip1 != NULL ) {
				computeDihedralAngle(N_i,CA_i,CO_i,N_ip1,psi);
			}
			if( CA_i != NULL && CO_i != NULL && N_ip1 != NULL 
				&& CA_ip1 != NULL ) {
				computeDihedralAngle(CA_i,CO_i,N_ip1,CA_ip1,omega);
			}
		}
		residue[i].setPhiPsiOmega(phi,psi,omega) ;
	}	
}


/* identifies hydrogen-bonded residues between backbone amine and carboxyl groups */
void CHAIN_t::compute_NHCO_HBond() {
	if( nResidues ==0 ) return ;
	double vwN = 1.55 ; // van der Waal's radius of N
	double vwO = 1.52 ; // van der Waal's radius of C
	double DELTA = 0.75 ; // const
	double HBNDLN = vwN+vwO+DELTA ;

	double vec1[3], vec2[3], theta ;
	for( uint32_t i = 1 ; i < nResidues ; i++ ) {
		double *N_i = NULL , *O_i = NULL, *C_i = NULL ; // N and C==O backbone atoms of ith res 
		 N_i = residue[i].getCoords(" N ") ;
		 O_i = residue[i].getCoords(" O ") ;
		 C_i = residue[i].getCoords(" C ") ;

		if( N_i == NULL ||  O_i == NULL || C_i == NULL ) continue ;
	
		double minhblen_d = INFINITY ;
		double minhblen_a = INFINITY ;
		bool hbD_setflag = false ;
		bool hbA_setflag = false ;
		for( uint32_t j = i+3 ; j < nResidues ; j++ ) { 
			double *N_j = NULL , *O_j = NULL, *C_j = NULL ; // N and C==O backbne atoms of jth res
			 N_j = residue[j].getCoords(" N ") ;
			 O_j = residue[j].getCoords(" O ") ;
		 	 C_j = residue[j].getCoords(" C ") ;

			if( N_j == NULL ||  O_j == NULL || C_j == NULL ) continue ;

			//testing (N-H)_i ----> (O=C)_j
			// NH group of i is the donar and CO group of j is the acceptor
			double distNiCOj  = normAminusB(N_i,O_j);
			//find angle between (C_j-O_j) and (N_i-O_j)
			vec1[0] = C_j[0]-O_j[0] ;
			vec1[1] = C_j[1]-O_j[1] ;
			vec1[2] = C_j[2]-O_j[2] ;
			
			vec2[0] = N_i[0]-O_j[0] ;
			vec2[1] = N_i[1]-O_j[1] ;
			vec2[2] = N_i[2]-O_j[2] ;
			theta = computeAngle( vec1, vec2 );
			if( distNiCOj <= HBNDLN && theta > 110 ) { 
				if(hbD_setflag == false) {
					//cout << chainID ;
					//cout << residue[i].getResIdentifier() 
					//	<< " bond donar to "
					//	<< residue[j].getResIdentifier() 
					//	<< " " << distNiCOj << "|" << HBNDLN
					//      << "|" << theta 
					//	<< endl ;
					residue[i].setNHCO_HBondResIndx(j,0) ;
					hbD_setflag = true ;
					minhblen_d = distNiCOj ;
				}
				else {
					if(distNiCOj <  minhblen_d)  {
					//	cout << chainID ;
					//	cout << residue[i].getResIdentifier() 
					//		<< " bond donar to "
					//		<< residue[j].getResIdentifier() 
					//		<< " " << distNiCOj << "|" << HBNDLN
					//      << "|" << theta 
					//		<< endl ;
						residue[i].setNHCO_HBondResIndx(j,0) ;
						minhblen_d = distNiCOj ;
					}

				}
			}


			//testing (C=0)_i ----> (H-N)_i
			// CO group of i is the acceptor and NH group of i is the donar
			double distCOiNj  = normAminusB(O_i,N_j);
			//find angle between (C_i-O_i) and (N_j-O_i)
			vec1[0] = C_i[0]-O_i[0] ;
			vec1[1] = C_i[1]-O_i[1] ;
			vec1[2] = C_i[2]-O_i[2] ;
			
			vec2[0] = N_j[0]-O_i[0] ;
			vec2[1] = N_j[1]-O_i[1] ;
			vec2[2] = N_j[2]-O_i[2] ;
			theta = computeAngle( vec1, vec2 );
			if( distCOiNj <= HBNDLN && theta > 110 ) { 
				if(hbA_setflag == false) {
					//cout << chainID ;
					//cout << residue[i].getResIdentifier() 
					//	<< " bond acceptor to "
					//	<< residue[j].getResIdentifier() 
					//	<< " " << distCOiNj << "|" << HBNDLN
					//      << "|" << theta 
					//	<< endl ;
					residue[i].setNHCO_HBondResIndx(j,1) ;
					hbA_setflag = true ;
					minhblen_a = distCOiNj ;
				}
				else {
					if(distCOiNj <  minhblen_a)  {
						//cout << chainID ;
						//cout << residue[i].getResIdentifier() 
						//<< " bond acceptor to "
						//<< residue[j].getResIdentifier() 
						//<< " " << distCOiNj << "|" << HBNDLN
					        //<< "|" << theta 
						//<< endl ;
						residue[i].setNHCO_HBondResIndx(j,1) ;
						minhblen_a = distNiCOj ;
					}

				}
			}

		}
	}	
}

double* CHAIN_t::getPhiPsiOmega( const uint32_t resIndx ) {
	assert( resIndx < nResidues ) ;
	return residue[resIndx].getPhiPsiOmega() ;
}

uint32_t * CHAIN_t::getNHCO_HBondResIndx( const uint32_t resIndx ) {
	assert( resIndx < nResidues ) ;
	return residue[resIndx].getNHCO_HBondResIndx() ;
}


void CHAIN_t::process( const char line[] ) {
	static char recordstr[7] ;
	static char resseqnumstr[6] ; // stores resSeq+Icode part of line
	static char aa_3letter[4] ;
	static map<string,char>::iterator aacitr;

	strncpy( recordstr, line, 6 ) ; recordstr[6] = '\0' ;
	
	aa_3letter[0] = toupper(line[17]) ;
	aa_3letter[1] = toupper(line[18]) ;
	aa_3letter[2] = toupper(line[19]) ;
	aa_3letter[3] = '\0' ;

	// parse only when the record type corresponds to a residue.
	// All residue entries start with an ATOM record. However in
	// some cases the residue MSE (Selenometheonine) starts with
	// HETATOM record type.
	if( strcmp( recordstr, "ATOM  ") == 0 
	    || (strcmp( recordstr, "HETATM") == 0
		&& strcmp( aa_3letter, "MSE")  == 0 )) {
		//copy resSeq+Icode information.
		resseqnumstr[0] = line[22] ;
		resseqnumstr[1] = line[23] ;
		resseqnumstr[2] = line[24] ;
		resseqnumstr[3] = line[25] ;
		resseqnumstr[4] = line[26] ;
		resseqnumstr[5] = '\0'     ;

		string tmp_aastr = aa_3letter ;
		char code ;
		
			
		aacitr = aacode.find(tmp_aastr) ;

		if( aacitr == aacode.end() ) {
      cerr << "STRANGE!!! " ;
			cerr << tmp_aastr << endl ;
		}
		// if the amino acid types is not among the list
		// declared in the constructor, then abort.
		assert( aacitr != aacode.end() ) ;
		
		//code is the single letter code of the aa type.
		code = aacitr->second ;

		int currResIndx = nResidues-1 ;

		//if new residue, create new residue object.
		if( isResidueIndexed( resseqnumstr ) == false ) {
			//insert into map, resseqnumstr-->vectorindx
			resNum_indx.insert( 
				pair<string,int>(resseqnumstr,nResidues) 
			);
			nResidues++ ;
			//append a new chain object
			RESIDUE_t newres ;
			residue.push_back( newres ) ;

			//assign aa_3letter name to the member in the object
			currResIndx = nResidues-1 ;
			residue[currResIndx].aaAssign( aa_3letter, code ) ;
		}
		//process the line into its corresponding residue.
		residue[currResIndx].process( line ) ;
	}
}

void CHAIN_t::setChainID( const char c )  {
	chainID = c ;
}

void CHAIN_t::print() {
	cout << endl ;
	cout << "\t\t" << "residue indx: \n" ;
	map<string,int>::iterator it ;

	for ( it=resNum_indx.begin() ; it != resNum_indx.end(); it++ )
		cout  << "\t\t" << (*it).first << " => " << (*it).second << endl;
	cout << endl ;

	cout  << "\t\t" << "nResidues = " << nResidues << endl ;
	for( unsigned int i = 0 ; i < nResidues ; i++ ) {
		cout  << "\t\t" << "RESIDUE " << i << endl ;
		residue[i].print() ;
	}
	//print sses and their types
	cout << "\t\t" << "nSSEs = " << nSSEs << endl ;
	for( uint32_t i = 0 ; i < nSSEs ; i++ ) {
		cout << "\t\t" << "Residues " << sse[i][0] ;
		cout << getResIdentifier(sse[i][0]) ;
		cout << " -- " ;
		cout << sse[i][1] ;
		cout << getResIdentifier(sse[i][1]) ;
		cout << " " << sseType[i] << endl ;
	}	

}

void CHAIN_t::printPDB( ofstream &out) {
	for( unsigned int i = 0 ; i < nResidues ; i++ ) {
		residue[i].printPDB(out) ;
	}
}

void CHAIN_t::printPDB( ofstream &out, Superpose3DClass &supobj) {
	for( unsigned int i = 0 ; i < nResidues ; i++ ) {
		residue[i].printPDB(out,supobj) ;
	}
}


// SSE related function defs
const uint32_t CHAIN_t::getnSSEs() const {
	if( sseFlag == false )  return 0 ;
	else return nSSEs ;
}

const vector<vector<uint32_t> >& CHAIN_t::getSSEBoundaries() const {
	//assert( sseFlag ) ;
	return sse ;
}

const vector<char>& CHAIN_t::getSSETypes() const {
	assert( sseFlag ) ;
	return sseType ;
}

char CHAIN_t::getSSType( const uint32_t resID ) const {
	return residue[resID].getSSType() ;
}
		

void CHAIN_t::setSSType( const uint32_t resID, const char c ) {
	residue[resID].setSSType(c) ;
}

/*appends both the SSE boundaries and type to respective vectors*/
void CHAIN_t::setSSE( const uint32_t s, const uint32_t e, const char c ) {
	assert( s < getnResidues() ) ;
	assert( e < getnResidues() ) ;

	if( sseFlag == false ) sseFlag = true ;
	vector<uint32_t> startAndEnd ;
	startAndEnd.push_back(s) ;
	startAndEnd.push_back(e) ;
	sse.push_back(startAndEnd) ;
	nSSEs++ ;
	sseType.push_back(c) ;
}

void CHAIN_t::sortSSE() {
	uint32_t n = sse.size() ;	
	if( n < 2 ) return ;

	assert( sseFlag ) ;
	uint32_t tmp[2] ;
	char tmp_c ;
	//bubble sort based on start points of sse
	uint32_t nprime ;
	do {
		nprime = 0 ;
		for( uint32_t i = 0 ; i < n-1 ; i++) {
			if( sse[i][0] > sse[i+1][0]) {
				//swap sse[i] and sse[i+1]
				tmp[0] = sse[i][0] ;
				tmp[1] = sse[i][1] ;
				sse[i][0] = sse[i+1][0] ;
				sse[i][1] = sse[i+1][1] ;
				sse[i+1][0] = tmp[0] ;
				sse[i+1][1] = tmp[1] ;
				//swap sseType[i] sseType[i+1]
				tmp_c = sseType[i] ;
				sseType[i] = sseType[i+1] ;
				sseType[i+1] = tmp_c ;
				nprime = i+1 ;
			}
		}
    		n = nprime ;
	} while(n > 1) ;
}

void CHAIN_t::genSSEPymolCmds(ofstream &of) {
	map<char,int>::iterator itr ;
	for( uint32_t i = 0 ; i < nSSEs ; i++ ) {
		of << "alter " 
			<< chainID
			<< "/"
			<< getResSeqNum(sse[i][0]) 
			<< "-"
			<< getResSeqNum(sse[i][1])
			<< "/, ss=" ;
		if( sseType[i] == 'G'
			|| sseType[i] == 'H'
			|| sseType[i] == 'I') {
			of << "\'H\'"<< endl ;
		}
		else if( sseType[i] == 'E'
			|| sseType[i] == 'L') {
			of << "\'S\'"<< endl ;
		}
		//color
		//first select
		stringstream tmps;
	       tmps << "sele" << chainID << i  ;
		//prepare selection name
		of << "select " 
			<< tmps.str() << ", "
			<< chainID << "/"
			<< getResSeqNum(sse[i][0]) 
			<< "-"
			<< getResSeqNum(sse[i][1])
			<< "/" << endl ;
		//color the above selection 
		//choose a random color
		int ci = rand()%nPymolColors ;
		of << "color " 
			<< pymolColor[ci]
			<< ", "
			<< tmps.str() 
			<< endl ;
	}	
}

void CHAIN_t::genSegmentPymolCmds(
  vector<vector<uint32_t> > &seg, 
  vector<char> &labels, 
  vector<string> &colors,
  ofstream &of
) {

	uint32_t n = seg.size() ;
	/*color each segment*/
        uint32_t cntrseg = 0 ;
	for( uint32_t i = 0 ; i < n ; i++ ) {
                //increment cntrseg if the segment is new (as against continuing
                //a previous one.
                if (seg[i][3] == 0) cntrseg++ ;
		stringstream locusstr;
		locusstr << chainID
			<< "/"
			<< getResSeqNum(seg[i][0]) 
			<< "-"
			<< getResSeqNum(seg[i][1])
			<< "/" ;
	
		stringstream selename;
		selename << cntrseg << "_" << labels[i] << "_" << chainID
			<< "_"
			<< getResSeqNum(seg[i][0]) 
			<< "-"
			<< getResSeqNum(seg[i][1]);

		//first select
		of << "select "  << selename.str() <<  ", "
				<< locusstr.str() << endl ;
		//then, color the above selection 
		of << "color " << colors[i] << ", "
				<< selename.str() << ", "
				<< endl ;
 
	}

	/* alter ss def based on def at residue level and 
	   NOT on the supplied Segment level labels */
	for( uint32_t i = 0 ; i < nResidues ; i++ ) {
		char sstype = residue[i].getSSType() ;
		//change 'C' defs to 'L'
		if(  sstype == 'C' || sstype == 'T' || sstype == '3' || sstype == '4' || sstype == '5' ) continue ;
		//change 'E' defs to 'S'
		else if(  sstype == 'E' ) sstype = 'S' ;
                //change to 'H' for all other cases
		else sstype = 'H' ;
		stringstream locusstr;
		locusstr << chainID
			<< "/"
			<< getResSeqNum(i) 
			<< "/" ;
		of << "alter " << locusstr.str() 
			<< ", ss=\'"<< sstype 
			<<"\'" << endl ; 
	}	
/*
	uint32_t n = seg.size() ;
	for( uint32_t i = 0 ; i < n ; i++ ) {
		stringstream locusstr;
		locusstr << chainID
			<< "/"
			<< getResSeqNum(seg[i][0]) 
			<< "-"
			<< getResSeqNum(seg[i][1])
			<< "/" ;

		stringstream selename;
		selename << i << "_" << labels[i] << "_" << chainID
			<< "_"
			<< getResSeqNum(seg[i][0]) 
			<< "-"
			<< getResSeqNum(seg[i][1]);

		if( labels[i] == 'A' || labels[i] == '3' 
			|| labels[i] == 'P' ) {
			of << "alter " << locusstr.str() 
				<< ", ss=\'H\'" << endl ; 

			//color
			//first select
			of << "select "  
				<< selename.str() <<  ", "
				<< locusstr.str() << endl ;
			//color the above selection 
			//choose a random color
			int ci = rand()%nPymolColors ;
			of << "color " 
				<< pymolColor[ci]
				<< ", "
				<< selename.str() << ", "
				<< endl ;
		}
		else if( labels[i] == 'E') {
			of << "alter " << locusstr.str() 
				<< ", ss=\'S\'" << endl ; 

			//color
			//first select
			of << "select "  
				<< selename.str() << ", "
				<< locusstr.str() << endl ;
			//color the above selection 
			//choose a random color
			int ci = rand()%nPymolColors ;
			of << "color " 
				<< pymolColor[ci]
				<< ", "
				<< selename.str() << ", "
				<< endl ;
		}
		else if( labels[i] == 'O' && 
			(seg[i][1]-seg[i][0]+1>=5) 
		) {
			of << "alter " << locusstr.str() 
				<< ", ss=\'L\'" << endl ; 

			//color
			//first select
			of << "select "  
				<< selename.str() <<  ", "
				<< locusstr.str() << endl ;
			//color the above selection 
			//choose a random color
			int ci = rand()%nPymolColors ;
			of << "color " 
				<< pymolColor[ci]
				<< ", "
				<< selename.str() << ", "
				<< endl ;
		}

	}
*/
}


void CHAIN_t::genSparseRepPDB(ofstream &of) {
	static string hetatm_prefix1 = "HETATM 2195  P   PO4   126    " ;
	static string hetatm_prefix2 = "HETATM 2196  O   PO4   126    " ;
	static string hetatm_suffix = "  0.50 59.85           P" ;

	assert(centersOfMass.size() == sse.size() ) ;
	double t ;		
	uint32_t n = centersOfMass.size() ;
	uint32_t outline_cntr = 0 ;
	char tstr[1000] ;
	for( uint32_t i = 0 ; i < n ; i++ ) {

		//first center of mass
		//1-6
		of << "HETATM" ;
		//7-11
		sprintf( tstr, "%5u", outline_cntr ) ;
		of << tstr ;
		//12-17
		of << "  P   " ;
		//18-20
		of << "P04" ;
		//21-22 blank
		of << "  " ;
		//23-26
		sprintf( tstr, "%4u", outline_cntr++ ) ;
		of << tstr ;
		//27-30 blank
		of << "    " ;

		t = centersOfMass[i][0] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;

		t = centersOfMass[i][1] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		
		t = centersOfMass[i][2] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		of << hetatm_suffix << endl ;

		/*
		//next another coordinate which is 
		//10A away from CM in the direction of SSE axis
		//1-6
		of << "HETATM" ;
		//7-11
		sprintf( tstr, "%5u", outline_cntr ) ;
		of << tstr ;
		//12-17
		of << "  O   " ;
		//18-20
		of << "OXY" ;
		//21-22 blank
		of << "  " ;
		//23-26
		sprintf( tstr, "%4u", outline_cntr++ ) ;
		of << tstr ;
		//27-30 blank
		of << "    " ;


		t = centersOfMass[i][0]+10*axes[i][0] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;

		t = centersOfMass[i][1]+10*axes[i][1];
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		
		t = centersOfMass[i][2]+10*axes[i][2]  ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		of << hetatm_suffix << endl ;
		*/


		//put atom at startpoint of the SSE.
		//1-6
		of << "HETATM" ;
		//7-11
		sprintf( tstr, "%5u", outline_cntr ) ;
		of << tstr ;
		//12-17
		of << "  O   " ;
		//18-20
		of << "OXY" ;
		//21-22 blank
		of << "  " ;
		//23-26
		sprintf( tstr, "%4u", outline_cntr++ ) ;
		of << tstr ;
		//27-30 blank
		of << "    " ;


		t = startPoints[i][0] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;

		t = startPoints[i][1] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		
		t = startPoints[i][2] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		of << hetatm_suffix << endl ;


		//put atom at endpoint of the SSE.
		//1-6
		of << "HETATM" ;
		//7-11
		sprintf( tstr, "%5u", outline_cntr ) ;
		of << tstr ;
		//12-17
		of << "  O   " ;
		//18-20
		of << "OXY" ;
		//21-22 blank
		of << "  " ;
		//23-26
		sprintf( tstr, "%4u", outline_cntr++ ) ;
		of << tstr ;
		//27-30 blank
		of << "    " ;


		t = endPoints[i][0] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;

		t = endPoints[i][1] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		
		t = endPoints[i][2] ;
		sprintf( tstr, "%8.3f", t ) ;
		of << tstr ;
		of << hetatm_suffix << endl ;
	}
	of.close() ;
}

bool CHAIN_t::computeAxisAndCM( const uint32_t sseindx) {
	static double axis[3] ;
	const uint32_t sp = sse[sseindx][0] ;
	const uint32_t ep = sse[sseindx][1] ;

	assert( sp < getnResidues() ) ;
	assert( ep < getnResidues() ) ;
	assert( ep-sp+1 >= 4 ) ;

	//extract CA coords from this region
	vector<vector<double> > CA ;
	for( uint32_t i = sp ;  i <= ep ;  i++ ) {
		double *tarr ;
		tarr = getCACoords(i) ;
		//mandatory check
		assert(  tarr != NULL ) ;

		vector<double> tvec ;
		tvec.push_back(tarr[0]) ;
		tvec.push_back(tarr[1]) ;
		tvec.push_back(tarr[2]) ;

		CA.push_back(tvec) ;
	}
	// Let 
	// CA[i-1], CA[i], and CA[i+1] be 3 position vectors. 
	// Compute:
	// v1 = CA[i-1] - CA[i] and 
	// v2 = CA[i+1] - CA[i]
	// Then compute midpoint as (v1+v2/2)
	vector<vector<double> > midpoints ;
	uint32_t n = CA.size() ;
	double v1[3] ;
	double v2[3] ;

	for( uint32_t i = 1 ; i < n-1 ;  i++ ) {
		v1[0] = CA[i-1][0]-CA[i][0] ;
		v1[1] = CA[i-1][1]-CA[i][1] ;
		v1[2] = CA[i-1][2]-CA[i][2] ;

		v2[0] = CA[i+1][0]-CA[i][0] ;
		v2[1] = CA[i+1][1]-CA[i][1] ;
		v2[2] = CA[i+1][2]-CA[i][2] ;

		vector<double> tmid;
		//find tmid = (v1+v2)/2
		tmid.push_back( (v1[0]+v2[0])/2 ) ;
		tmid.push_back( (v1[1]+v2[1])/2 ) ;
		tmid.push_back( (v1[2]+v2[2])/2 ) ;
		//convert tmid to position vector relative to CA[i]
		tmid[0] += CA[i][0] ;
		tmid[1] += CA[i][1] ;
		tmid[2] += CA[i][2] ;

		midpoints.push_back( tmid ) ;
	}

	//compute center of mass of all the midpoints
	vector<double> cm ;
	cm.push_back(0) ;
	cm.push_back(0) ;
	cm.push_back(0) ;

	n = midpoints.size() ;

	for( uint32_t i = 0 ; i < n ;  i++ ) {
		cm[0] += midpoints[i][0] ;
		cm[1] += midpoints[i][1] ;
		cm[2] += midpoints[i][2] ;
	}
	cm[0]/=n ; cm[1]/=n ; cm[2]/=n ;

	// copy local cm to main data structure.
	centersOfMass.push_back( cm ) ;

	/* Now build the matrix A as a nx3 matrix, where 
	 * the first column is  x_i - xbar
	 * the second column is y_i - ybar
	 * the third column is  z_i - zbar
	 */
	vector<vector<double> > A ;
	for( uint32_t i = 0 ; i < n ; i++ ) {
		vector<double> tmp ;
		tmp.push_back( midpoints[i][0] - cm[0] ) ;
		tmp.push_back( midpoints[i][1] - cm[1] ) ;
		tmp.push_back( midpoints[i][2] - cm[2] ) ;
		A.push_back( tmp ) ;
	}
	double w[3] ; // to store singular values in decreasing order
	double v[3][3] ; // store corresponding singular vectors

	//SV decompose matrix A
	int status = fsvd( A, w, v ) ;

	/* if method fails use the alternative method that 
	* defines the axis as the line joining the top
	* and bottom midpoints.
	*/
	if( status == 0 ) {
		double norm = normAminusB(midpoints[n-1],midpoints[0]);
		axis[0]=(midpoints[n-1][0]-midpoints[0][0])/norm;
		axis[1]=(midpoints[n-1][1]-midpoints[0][1])/norm;
		axis[2]=(midpoints[n-1][2]-midpoints[0][2])/norm;
	}
	/* else if the svd was successful*/
	else {

		// the direction cosines of the axis are given by the 
		// vector corresponding to the smallest singular value
		axis[0] = v[0][0] ;
		axis[1] = v[1][0] ;
		axis[2] = v[2][0] ;
		//cout << "norm : " << vectorNorm( axis ) << endl  ;

		// However, we have to check the direction of axis so 
		// that it is along N- to C- direction and not the other
		// way. 
		// Compute dot product with the approximate 
		// N- to C- midpoint  vector and check the sign

		double taxis[3] ;
		double norm = normAminusB( midpoints[n-1], midpoints[0] ) ;
		taxis[0] = (midpoints[n-1][0]-midpoints[0][0])/norm ;
		taxis[1] = (midpoints[n-1][1]-midpoints[0][1])/norm ;
		taxis[2] = (midpoints[n-1][2]-midpoints[0][2])/norm ;

		double dprod = -999.999 ;
		computeDotProduct(taxis, axis, dprod ) ;
		if( dprod < 0 ) {
			// this means that the smalled singular vecotr 
			// is pointing the other way from C- to N- 
			// direction. Reorient the axis to point in the 
			// N- to C- direction
			axis[0] *= -1 ;
			axis[1] *= -1 ;
			axis[2] *= -1 ;
		}
	}
	//prepare to return a std::vector form of the axis
	vector<double> axis_vec ; 
	axis_vec.push_back( axis[0] ) ;
	axis_vec.push_back( axis[1] ) ;
	axis_vec.push_back( axis[2] ) ;
	axes.push_back( axis_vec ) ;


	// This code has been added subsequently.
	// We will now compute the start and end points of the SSE.
	// This is done by projecting the start and end residue 
	// coordinates onto the axis.
	double pA[3], pB[3], V[3], pP[3] ;

	//project startpoint
	double *tarr = getCACoords(sp) ;
	pA[0] = tarr[0] ;
	pA[1] = tarr[1] ;
	pA[2] = tarr[2] ;

	pB[0] = cm[0] ; 
	pB[1] = cm[1] ; 
	pB[2] = cm[2] ; 

	V[0] = axis[0] ;
	V[1] = axis[1] ;
	V[2] = axis[2] ;
	projectPoint2Line(pA,pB,V,pP) ;
	vector<double> sp_vec ; 
	sp_vec.push_back( pP[0] ) ;
	sp_vec.push_back( pP[1] ) ;
	sp_vec.push_back( pP[2] ) ;
	startPoints.push_back( sp_vec ) ;

	tarr = getCACoords(ep) ;
	pA[0] = tarr[0] ;
	pA[1] = tarr[1] ;
	pA[2] = tarr[2] ;

	pB[0] = cm[0] ; 
	pB[1] = cm[1] ; 
	pB[2] = cm[2] ; 

	V[0] = axis[0] ;
	V[1] = axis[1] ;
	V[2] = axis[2] ;
	projectPoint2Line(pA,pB,V,pP) ;

	vector<double> ep_vec ; 
	ep_vec.push_back( pP[0] ) ;
	ep_vec.push_back( pP[1] ) ;
	ep_vec.push_back( pP[2] ) ;
	endPoints.push_back( ep_vec ) ;

	return true ;
}


void CHAIN_t::computeAxesAndCentersOfMass() {
	const uint32_t n = sse.size() ;	

	for( uint32_t i = 0 ; i < n ; i++ ) {
		bool status = computeAxisAndCM(i) ;
		assert( status ) ;
	}
	axesAndCMFlag = true ;
}

const vector<vector<double> >& CHAIN_t::getSSEAxes() const {
	assert( axesAndCMFlag ) ;
	return axes ;
}


const vector<vector<double> >& CHAIN_t::getSSECentersOfMass() const {
	assert( axesAndCMFlag ) ;
	return centersOfMass ;
}

const vector<vector<double> >& CHAIN_t::getSSEStartPoints() const {
	assert( axesAndCMFlag ) ;
	return startPoints ;
}

const vector<vector<double> >& CHAIN_t::getSSEEndPoints() const {
	assert( axesAndCMFlag ) ;
	return endPoints ;
}


/* These are the mapping used for various conversions between
 * orientation angles, double-quandrant encoding and bytecodes
 * ---------------------------------------------------------
 * Omega	Tableau Encoding	ByteCode
 * [0,45) 		PD 		0
 * [45,90) 		RD 		1
 * [90,135) 		RT 		2
 * [135,180) 		OT 		3
 *
 * [-180,-135) 		OS 		4
 * [-135,-90) 		LS 		5
 * [-90,-45) 		LE 		6
 * [-45,0) 		PE 		7
 * ---------------------------------------------------------
 * */
/* tableau related function*/
void CHAIN_t::omegaToBytecode( const double omega, uint8_t &code) {
	//special case, as I am not sure if my geometry3D routine
	//gives 180 or -180 
	if( omega == 180 ) {
		code = 4 ;
		return ;
	}

	if( omega >= 0 ) {
		if( omega >= 0 && omega < 45 ) {
			code = 0 ;
		}
		else if( omega >= 45 && omega < 90 ) {
			code = 1 ;
		}
		else if( omega >= 90 && omega < 135 ) {
			code = 2 ;
		}
		else if( omega >= 135 && omega < 180 ) {
			code = 3 ;
		}
	}
	else {
		if( omega >= -180 && omega < -135 ) {
			code = 4 ;
		}
		else if( omega >= -135 && omega < -90 ) {
			code = 5 ;
		}
		else if( omega >= -90 && omega < -45 ) {
			code = 6 ;
		}
		else if( omega >= -45 && omega < 0 ) {
			code = 7 ;
		}
	}
}

void CHAIN_t::omegaToEncoding(  const double omega, char encode[2] ) {
	//special case, as I am not sure if my geometry3D routine
	//gives 180 or -180 
	if( omega == 180 ) {
		encode[0] = 'O' ;
		encode[1] = 'S' ;
		return ;
	}

	if( omega >= 0 ) {
		if( omega >= 0 && omega < 45 ) {
			encode[0] = 'P' ;
			encode[1] = 'D' ;
		}
		else if( omega >= 45 && omega < 90 ) {
			encode[0] = 'R' ;
			encode[1] = 'D' ;
		}
		else if( omega >= 90 && omega < 135 ) {
			encode[0] = 'R' ;
			encode[1] = 'T' ;
		}
		else if( omega >= 135 && omega < 180 ) {
			encode[0] = 'O' ;
			encode[1] = 'T' ;
		}
	}
	else {
		if( omega >= -180 && omega < -135 ) {
			encode[0] = 'O' ;
			encode[1] = 'S' ;
		}
		else if( omega >= -135 && omega < -90 ) {
			encode[0] = 'L' ;
			encode[1] = 'S' ;
		}
		else if( omega >= -90 && omega < -45 ) {
			encode[0] = 'L' ;
			encode[1] = 'E' ;
		}
		else if( omega >= -45 && omega < 0 ) {
			encode[0] = 'P' ;
			encode[1] = 'E' ;
		}
	}
}

void CHAIN_t::byteCodeToEncoding( const uint8_t code, char encode[2] ) {
	switch( code ) {
		case 0:
			encode[0] = 'P' ;
			encode[1] = 'D' ;
			break;
		case 1:
			encode[0] = 'R' ;
			encode[1] = 'D' ;
			break;
		case 2:
			encode[0] = 'R' ;
			encode[1] = 'T' ;
			break;
		case 3:
			encode[0] = 'O' ;
			encode[1] = 'T' ;
			break;
		case 4:
			encode[0] = 'O' ;
			encode[1] = 'S' ;
			break;
		case 5:
			encode[0] = 'L' ;
			encode[1] = 'S' ;
			break;
		case 6:
			encode[0] = 'L' ;
			encode[1] = 'E' ;
			break;
		case 7:
			encode[0] = 'P' ;
			encode[1] = 'E' ;
			break;
	}
}

void CHAIN_t::computeTableau() {
	if(sse.size() <= 1) return ;
	assert( sseFlag == true ) ;	
	assert( axesAndCMFlag == true ) ;	

	const uint32_t nSSEs = sse.size() ;
	// local omega and tableau
	vector<vector<double> > omega(nSSEs, 
			vector<double>(nSSEs,-999.99) ) ;
	vector<vector<uint8_t> > tableau(nSSEs,
			vector<uint8_t>(nSSEs,255) );

	for( uint32_t i = 0 ; i < nSSEs ; i++ ) {
		for( uint32_t j = i+1 ; j < nSSEs ; j++ ) {
			double cmA[3], axisA[3],cmB[3],axisB[3] ;
			cmA[0] = centersOfMass[i][0] ;
			cmA[1] = centersOfMass[i][1] ;
			cmA[2] = centersOfMass[i][2] ;

			axisA[0] = axes[i][0] ;
			axisA[1] = axes[i][1] ;
			axisA[2] = axes[i][2] ;

			cmB[0] = centersOfMass[j][0] ;
			cmB[1] = centersOfMass[j][1] ;
			cmB[2] = centersOfMass[j][2] ;

			axisB[0] = axes[j][0] ;
			axisB[1] = axes[j][1] ;
			axisB[2] = axes[j][2] ;
			omega[i][j] = getOrientationAngle(cmA,axisA,cmB,axisB);
			omega[j][i] = omega[i][j] ;

			uint8_t code ;
			omegaToBytecode( omega[i][j], code ) ; 
			tableau[i][j] = code ; 
			tableau[j][i] = tableau[i][j] ; 
		}
	}
	tableauFlag = true ;
	this->omega = omega ;
	this->tableau = tableau ;
}

const vector<vector<double> >& CHAIN_t::getTableau_orientations() const {
	assert( tableauFlag) ;
	return omega ;
}

const vector<vector<uint8_t> >& CHAIN_t::getTableau_codes() const {
	assert( tableauFlag) ;
	return tableau ;
}

void CHAIN_t::printTableauInfo(ofstream &out) {
	const uint32_t n = sse.size() ;
	if(n < 2 ) return ;
	out << "list of SSEs\n" ;
	for( uint32_t i = 0 ; i < n ; i ++ ) {
		out << i << " - " 
			<< sse[i][0] <<  getResIdentifier(sse[i][0])
			<< " to "
			<< sse[i][1] <<  getResIdentifier(sse[i][1])
			<< " is a "
			<< sseType[i] << endl ;
	}

	out << "omega matrix:\n" ;
	for( uint32_t i = 0 ; i < n ; i ++ ) {
		for( uint32_t j = 0 ; j < n; j++ ) {
			out << fixed ;
			out << setprecision(1) ;
			if(i != j ){
				out << setw(8) ;
				out << omega[i][j] ;
			}
			else {
				out << "     *" ;
				out << sseType[j] ;
				out << "*" ;

			}
		}
		out << endl ;
	}
	out << endl ;
	out << "tableau matrix:\n" ;
	for( uint32_t i = 0 ; i < n ; i ++ ) {
		for( uint32_t j = 0 ; j < n; j++ ) {
			out << fixed ;
			out << setprecision(1) ;
			char encode[2] ;
			byteCodeToEncoding( tableau[i][j], encode ) ;
			if(i != j )
				out << " " << encode[0] << encode[1] << " " ; 
			else
				out << " ** " ; 
		}
		out << endl ;
	}
	out << endl ;
}



/*#### MODEL_t CLASS ###*/
MODEL_t::MODEL_t() {
		  nChains = modelID = 0 ;
}

const uint32_t MODEL_t::getnChains() const {
	return nChains ;
}

const uint32_t MODEL_t::getnResidues( const uint32_t chainIndx ) const {
	assert( chainIndx < nChains ) ; 
	return chain[chainIndx].getnResidues() ; 
}

const uint32_t MODEL_t::getnAtoms( 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) const 
{
	assert( chainIndx < nChains ) ; 
	return chain[chainIndx].getnAtoms( resIndx ) ; 
}

char MODEL_t::getResCode( 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) const 
{
	assert( chainIndx < nChains ) ; 
	return chain[chainIndx].getResCode( resIndx ) ; 
}

string MODEL_t::getResIdentifier( 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) const 
{
	assert( chainIndx < nChains ) ; 
	return chain[chainIndx].getResIdentifier( resIndx ) ; 
}

int32_t MODEL_t::getResSeqNum( 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) const 
{
	assert( chainIndx < nChains ) ; 
	return chain[chainIndx].getResSeqNum(resIndx ) ; 
}


const ATOM_t& MODEL_t::getAtomRecord( 
  const uint32_t chainIndx, 
  const uint32_t resIndx, 
  const char pattern[] 
) {
	assert( chainIndx < nChains ) ; 
	return chain[chainIndx].getAtomRecord( resIndx, pattern ) ; 
}

double* MODEL_t::getCACoords( 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) {
	assert( chainIndx < nChains ) ; 
	return chain[chainIndx].getCACoords( resIndx ) ; 
}

char MODEL_t::getChainID( const uint32_t chainIndx ) {
	return chain[chainIndx].getChainID() ;
}

string MODEL_t::getAminoAcidSequenceOfChain(const uint32_t chainIndx) {
  return chain[chainIndx].getAminoAcidSequenceOfChain();
}

vector<string> MODEL_t::getResIdentifiersOfChain(const uint32_t chainIndx) {
  return chain[chainIndx].getResIdentifiersOfChain();
}

vector<vector<double> > MODEL_t::getCACoordsOfChain(const uint32_t chainIndx) {
  return chain[chainIndx].getCACoordsOfChain();
}

void MODEL_t::getCACoordsOfChain(
  const uint32_t chainIndx,
  vector<vector<double> > &CA,
  string &aaseq,
  vector<string> &resIDStrings,
  vector<uint32_t> &resObjectIDs
) {
  return chain[chainIndx].getCACoordsOfChain(CA,aaseq,resIDStrings,resObjectIDs);
}

void MODEL_t::getCACoordsOfChain(
  const char chainSymb,
  vector<vector<double> > &CA,
  string &aaseq,
  vector<string> &resIDStrings,
  vector<uint32_t> &resObjectIDs
) {
  uint32_t chainIndx;
  assert(isChainIndexed(chainSymb,chainIndx) == true);
  return chain[chainIndx].getCACoordsOfChain(CA,aaseq,resIDStrings,resObjectIDs);
}



void MODEL_t::computePhiPsiOmega() {
	for( uint32_t i = 0 ; i < nChains ; i++ ) {
		chain[i].computePhiPsiOmega() ;
	}
}

void MODEL_t::compute_NHCO_HBond() {
	for( uint32_t i = 0 ; i < nChains ; i++ ) {
		chain[i].compute_NHCO_HBond() ;
	}
}


bool MODEL_t::isChainIndexed(const char chainID, uint32_t &chainIndx) {
		static map<char,int>::iterator itr ;
		itr = chainID_indx.find( chainID );
      chainIndx = itr->second;

		if( itr == chainID_indx.end() ) {
         return false ;
		}
		else return true ;
}

void MODEL_t::process( const char line[] ) {
	static char recordstr[7] ;
	static char chainID ;
   uint32_t chainIndx;

	strncpy( recordstr, line, 6 ) ; recordstr[6] = '\0' ;
	if( strcmp( recordstr, "ATOM  ") == 0 
		|| strcmp( recordstr, "HETATM") == 0 
	) {
		chainID = line[21] ;
		// If this chainID is new, then append a new
		// chain object into the model.
		if( isChainIndexed(chainID,chainIndx) == false ) {
			//insert into map, chainID-->vectorindx
			chainID_indx.insert(pair<char,int>(chainID,nChains));
			nChains++ ;
			//append a new chain object
			CHAIN_t newchain ;
			newchain.setChainID( chainID ) ;
			chain.push_back( newchain ) ;
		}
		//process the line into its corresponding chain.
		int currChainIndx = nChains-1 ;
		chain[currChainIndx].process( line ) ;
	}
}


/*This returns a pointer to a double[3] containing dihedrals Phi, Psi, Omega */
double* MODEL_t::getPhiPsiOmega(
  const uint32_t chainIndx, 
  const uint32_t resIndx 
)  
{
	assert( chainIndx < nChains ) ;
	return chain[chainIndx].getPhiPsiOmega( resIndx ) ; 
}

uint32_t* MODEL_t::getNHCO_HBondResIndx(
  const uint32_t chainIndx, 
  const uint32_t resIndx 
)  
{
	assert( chainIndx < nChains ) ;
	return chain[chainIndx].getNHCO_HBondResIndx( resIndx ) ; 
}


void MODEL_t::print() {
	cout << endl ;
	cout << "\t" << "Chain labels: \n" ;
	map<char,int>::iterator it ;

	for ( it=chainID_indx.begin() ; it != chainID_indx.end(); it++ )
		cout << "\t"<< (*it).first << " => " << (*it).second << endl;
	cout << endl ;

	cout << "\t"<< "nChains = " << nChains << endl ;
	for( unsigned int i = 0 ; i < nChains ; i++ ) {
		cout<< "\t" << "CHAIN " << i << endl ;
		chain[i].print() ;
	}
}

void MODEL_t::printPDB(ofstream &out) {
	for( unsigned int i = 0 ; i < nChains ; i++ ) {
		chain[i].printPDB(out) ;
    out << "TER" << endl;
	}
}

void MODEL_t::printPDB(ofstream &out, Superpose3DClass &supobj) {
	for( unsigned int i = 0 ; i < nChains ; i++ ) {
		chain[i].printPDB(out,supobj) ;
    out << "TER" << endl;
	}
}



// secondary structure related functions
const uint32_t MODEL_t::getnSSEs( const uint32_t chainID) {
	return chain[chainID].getnSSEs() ;
}

const vector<vector<uint32_t> >& MODEL_t::getSSEBoundaries(
  const uint32_t chainID
) const {
	return chain[chainID].getSSEBoundaries() ;
}

const vector<char>& MODEL_t::getSSETypes(
  const uint32_t chainID
) const {
	return chain[chainID].getSSETypes() ;
}

char MODEL_t::getSSType(
  const uint32_t chainID,
  const uint32_t resID
) const{
	 
	return	chain[chainID].getSSType( resID);

}
		

void MODEL_t::setSSType( 
  const uint32_t chainID, 
  const uint32_t resID, 
  const char c 
) {
	chain[chainID].setSSType(resID,c) ;
}



/*appends both the SSE boundary and type*/
void MODEL_t::setSSE( 
  const uint32_t chainID,
  const uint32_t s, // start point
  const uint32_t e, // end point
  const char c // SSE type char code 
) {
	chain[chainID].setSSE(s,e,c) ;
}

void MODEL_t::sortSSE() {
	const uint32_t nC = getnChains() ;	
	for( uint32_t c = 0 ; c < nC ; c++ ) {
		chain[c].sortSSE() ;
	}
}

void MODEL_t::genSSEPymolCmds(ofstream &of ) {
	for( uint32_t i = 0 ; i < nChains ; i ++  ) {
		chain[i].genSSEPymolCmds(of) ;
	}

}

void MODEL_t::genSegmentPymolCmds(
  vector<vector<uint32_t> > &seg, 
  vector<char> &labels, 
  vector<string> &colors,
  uint32_t chainIndx,
  ofstream &of ) {
	assert( chainIndx < nChains ) ; 
	chain[chainIndx].genSegmentPymolCmds(seg,labels,colors, of) ;
}

void MODEL_t::genSparseRepPDB(ofstream &of ) {
	for( uint32_t i = 0 ; i < nChains ; i ++  ) {
		chain[i].genSparseRepPDB(of) ;
	}
}

void MODEL_t::computeAxesAndCentersOfMass() {
	const uint32_t nC = getnChains() ;	
	for( uint32_t c = 0 ; c < nC ; c++ ) {
		chain[c].computeAxesAndCentersOfMass() ;
	}
}

const vector<vector<double> >& MODEL_t::getSSEAxes(
  const uint32_t chainID
) const {
	return chain[chainID].getSSEAxes() ;
}

const vector<vector<double> >& MODEL_t::getSSECentersOfMass(
  const uint32_t chainID
) const {
	return chain[chainID].getSSECentersOfMass() ;
}

const vector<vector<double> >& MODEL_t::getSSEStartPoints(
  const uint32_t chainID
) const {
	return chain[chainID].getSSEStartPoints() ;
}

const vector<vector<double> >& MODEL_t::getSSEEndPoints(
  const uint32_t chainID
) const {
	return chain[chainID].getSSEEndPoints() ;
}

/* tableau related function*/
void MODEL_t::computeTableau() {
	const uint32_t nC = getnChains() ;	
	for( uint32_t c = 0 ; c < nC ; c++ ) {
		chain[c].computeTableau() ;
	}
}

const vector<vector<double> >& MODEL_t::getTableau_orientations(
  const uint32_t chainID
) const {
	return chain[chainID].getTableau_orientations() ;
}

const vector<vector<uint8_t> >& MODEL_t::getTableau_codes(
  const uint32_t chainID
) const {
	return chain[chainID].getTableau_codes() ;
}

void MODEL_t::printTableauInfo(ofstream &out ) {
	for( uint32_t i = 0 ; i < nChains ; i ++  ) {
		out << "+++++ CHAIN " << i << " +++++" << endl ;
		chain[i].printTableauInfo(out) ;
	}
}




/*### PDB_t CLASS ###*/
PDB_t::PDB_t( const string fname ) {
	nModels = 0 ;
	if( parsePDB( fname ) == false ) {
		cerr << "::PDB_PARSE_ERROR:: The file: " << fname  
			<< " is either NOT a PDB format file, or "
			<< " does NOT conform to its basic requirements.\n";
		exit(1) ;
	}
	computePhiPsiOmega() ;
	this->fname = fname ;
}

const string PDB_t::getTitle() const {
	return title ;
}

const string PDB_t::getTitleSummary() const {
	return title.substr(10,50) ;
}

string PDB_t::getFileName() const {
	return fname ;
}

const uint32_t PDB_t::getnModels() const {
	return nModels ;
}

const uint32_t PDB_t::getnChains( const uint32_t modelIndx ) const {
	assert( modelIndx < nModels ) ;

	return model[modelIndx].getnChains() ;
}

char PDB_t::getResCode( 
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) const 
{
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getResCode( chainIndx, resIndx ) ; 
}

/* This returns a concetenation of resi name and number*/
string PDB_t::getResIdentifier( 
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) const 
{
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getResIdentifier( chainIndx, resIndx ) ; 
}

/*This returns the residue sequence number*/
int32_t PDB_t::getResSeqNum( 
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
) const 
{
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getResSeqNum( chainIndx, resIndx ) ; 
}


/*This returns a pointer to a double[3] containing dihedrals Phi, Psi, Omega */
double* PDB_t::getPhiPsiOmega(
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
)  
{
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getPhiPsiOmega( chainIndx, resIndx ) ; 
}

uint32_t* PDB_t::getNHCO_HBondResIndx(
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t resIndx 
)  
{
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getNHCO_HBondResIndx( chainIndx, resIndx ) ; 
}



const ATOM_t& PDB_t::getAtomRecord( 
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t resIndx,
  const char pattern[] )  
{
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getAtomRecord(chainIndx,resIndx, pattern ) ; 
}

double* PDB_t::getCACoords( 
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t resIndx
)  {
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getCACoords(chainIndx,resIndx ) ; 
}

const uint32_t PDB_t::getnResidues(
  const uint32_t modelIndx, 
  const uint32_t chainIndx
) const
{
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getnResidues( chainIndx) ;
}

const uint32_t PDB_t::getnAtoms(
  const uint32_t modelIndx, 
  const uint32_t chainIndx, 
  const uint32_t atomIndx
) const {
	assert( modelIndx < nModels ) ;
	return model[modelIndx].getnAtoms( chainIndx, atomIndx ) ;
}

char PDB_t::getChainID(
  const uint32_t modelIndx, 
  const uint32_t chainIndx
) {
	return model[modelIndx].getChainID(chainIndx) ;
}

string PDB_t::getAllChainIDs(const uint32_t modelIndx){
   string allchainids;
   for (uint32_t chainIndx = 0; chainIndx < getnChains(modelIndx); chainIndx++){
	   allchainids.push_back(model[modelIndx].getChainID(chainIndx));
   }
   return allchainids;
}


string PDB_t::getAminoAcidSequenceOfChain(
  const uint32_t modelIndx,
  const uint32_t chainIndx
) {
  return model[modelIndx].getAminoAcidSequenceOfChain(chainIndx);
}

vector<string> PDB_t::getResIdentifiersOfChain(
  const uint32_t modelIndx,
  const uint32_t chainIndx
) {
  return model[modelIndx].getResIdentifiersOfChain(chainIndx);
}

vector<vector<double> > PDB_t::getCACoordsOfChain(
  const uint32_t modelIndx,
  const uint32_t chainIndx
) {
  return model[modelIndx].getCACoordsOfChain(chainIndx);
}

void PDB_t::getCACoordsOfChain(
  const uint32_t modelIndx,
  const uint32_t chainIndx,
  vector<vector<double> > &CA,
  string &aaseq,
  vector<string> &resIDStrings,
  vector<uint32_t> &resObjectIDs
) {
  model[modelIndx].getCACoordsOfChain(chainIndx,CA,aaseq,resIDStrings,resObjectIDs);
}

void PDB_t::getCACoordsOfChain(
  const uint32_t modelIndx,
  const char chainSymb,
  vector<vector<double> > &CA,
  string &aaseq,
  vector<string> &resIDStrings,
  vector<uint32_t> &resObjectIDs
) {
  model[modelIndx].getCACoordsOfChain(chainSymb,CA,aaseq,resIDStrings,resObjectIDs);
}

void PDB_t::getCACoordsOfChains(
  const uint32_t modelIndx,
  const string chainIDs,
  vector<vector<vector<double> > > &CAs,
  vector<string> &aaseqs,
  vector<vector<string> >&resIDStrings,
  vector<vector<uint32_t> > &resObjectIndxs
) {
   assert(chainIDs.length() > 0);
   assert(chainIDs.length() == aaseqs.size());
   for (uint32_t i = 0; i < chainIDs.length(); i++) {
      model[modelIndx].getCACoordsOfChain(chainIDs[i],CAs[i],aaseqs[i],resIDStrings[i],resObjectIndxs[i]);
   }
}






void PDB_t::computePhiPsiOmega() {
	for( uint32_t i = 0 ; i < nModels ; i++ ) {
		model[i].computePhiPsiOmega() ;
	}
}

void PDB_t::compute_NHCO_HBond() {
	for( uint32_t i = 0 ; i < nModels ; i++ ) {
		model[i].compute_NHCO_HBond() ;
	}
}


bool PDB_t::parsePDB( const string fname ) {
	// open pdb
	ifstream infile( fname.c_str(), ios::in );
	if(!infile) {
		cerr << "::PDB_PARSE_ERROR:: Unable to find the file: " 
			<< fname << "\n" ;
		exit(1) ;
	}

	char line[1000] ;
	char recordstr[7] ;
	bool headerZoneFlag = true ;
	uint32_t linecnt = 0 ; // line number counter
	while( !infile.eof()  ) {
		infile.getline( line, 1000 ) ;
		linecnt++ ; 

		if( infile.eof() ) break ; 
		if( strcmp( line, "" ) == 0 ) continue ;
		
		/* check whether the line starts with a valid RECORD TYPE */
		strncpy( recordstr, line, 6 ) ; recordstr[6] = '\0' ;
		if( isValidRecord( recordstr) == false ) {
		  cerr << "::PDB_PARSE_ERROR:: The parser stumbled into an "
			 << "unrecognized record type: " << recordstr 
			 << "at Line " << linecnt << " in File: " << fname
			 << "\n" ;
			 return false ;
		}

		if( strcmp( recordstr, "TITLE ") == 0 ) {
			char tline[1000] ;
			strcpy( tline, line ) ;
			int i ;
			for( i = strlen(tline)-1 ; i >= 0; i-- ) {
				if( tline[i] !=' ' ) break ;
			}
			tline[i+1] = '\0' ;	
			string tstr = tline ;
			title.append(tstr.substr(10) ) ;
		}

		// When a new MODEL line is encountered,
		// create a new model in the pdb
		if( strcmp( recordstr, "MODEL ") == 0 ) {
			if( headerZoneFlag == true) headerZoneFlag = false ;
			//create a new model object
			nModels++ ;
			MODEL_t newmodel ;
			model.push_back( newmodel ) ;
			continue ;
		} 

		// Handling the case where a MODEL line was not encountered
		// before encountering a ATOM or HETATM line. This happens 
		// in PDB files without MODEL records -- usually X-ray 
		// determined data (as against NMR). Treat such files
		// as a pdb file containing 1 MODEL.
		if( headerZoneFlag == true 
			&& (strcmp( recordstr, "ATOM  ") == 0 
				|| strcmp( recordstr, "HETATM") == 0 
			)
		) {
			headerZoneFlag = false ;
			//create a new model object
			nModels++ ;
			MODEL_t newmodel ;
			model.push_back( newmodel ) ;
		}

		if( headerZoneFlag == true ) continue ;
		else {
		/* the line that is to be parsed belongs to a model that 
		 * was previously created. The index of the object is equal
		 * to nModels-1. Process this line as a record of that
		 * object. */
			int currModelIndx = nModels-1 ;
			model[currModelIndx].process( line ) ;
		}
	}
	infile.close() ;
	return true ;
}

void PDB_t::print() {
	cout << "nModels = " << nModels << endl ;
	for( unsigned int i = 0 ; i < nModels ; i++ ) {
		cout << "MODEL " << i << endl ;
		model[i].print() ;
	}
}

void PDB_t::printPDB(string fname) {
  ofstream out(fname.c_str(),ios::out);
  assert(out);
	for( unsigned int i = 0 ; i < nModels ; i++ ) {
    if (nModels > 1) out << "MODEL " << i+1 << endl;
		model[i].printPDB(out) ;
    if (nModels > 1) out << "ENDMDL" << endl;
	}
  out.close();
}

void PDB_t::printPDB(string fname, Superpose3DClass &supobj) {
  ofstream out(fname.c_str(),ios::out);
  assert(out);
	for( unsigned int i = 0 ; i < nModels ; i++ ) {
    if (nModels > 1) out << "MODEL " << i+1 << endl;
		model[i].printPDB(out,supobj) ;
    if (nModels > 1) out << "ENDMDL" << endl;
	}
  out.close();
}



bool PDB_t::isValidRecord( const char str[] ) {
	static int nRecords = 69 ;
	static char recordType[69][7] =  {
		"ATOM  ", "ANISOU", "HETATM", "HEADER", "REMARK", 
		"AUTHOR", "CAVEAT", "COMPND", "EXPDTA", "MDLTYP", 
		"KEYWDS", "OBSLTE", "SOURCE", "SPLIT ", "SPRSDE", 
		"TITLE ", "CISPEP", "CONECT", "DBREF ", "HELIX ", 
		"SHEET ", "HET   ", "LINK  ", "MODRES", "MTRIX1", 
		"MTRIX2", "MTRIX3", "REVDAT", "SEQADV", "SHEET" , 
		"SSBOND", "FORMUL", "HETNAM", "HETSYN", "SEQRES", 
		"SITE  ", "ENDMDL", "MODEL ", "TER   ", "JRNL  ", 
		"SPLIT ", "NUMMDL", "DBREF1", "DBREF2", "HET   ", 
		"CRYST1", "CRYST2", "CRYST3", "ORIGX1", "ORIGX2", 
		"ORIGX3", "SCALE1", "SCALE2", "SCALE3", "MASTER", 
		"END   ", "FTNOTE", "TURN  ", "TVECT ",
		//shorter version of record types. Illegal but
		//allows parsing of some format variants
		"END  ",  "END ",   "END", "TER  ", "TER ", 
		"TER",
		// other records that were encountered while 
		// processing protein data bank
		"HYDBND", "SLTBRG", "SIGATM", "SIGUIJ"

	} ;
	

	for( int i = 0 ; i < nRecords ; i++ ) {
		if( strcmp( str, recordType[i] ) == 0 ) return true ;
	}
	return false ;
}

// secondary structure related functions
const vector<vector<uint32_t> >& PDB_t::getSSEBoundaries(
  const uint32_t modelID, 
  const uint32_t chainID
) const {
	return model[modelID].getSSEBoundaries(chainID) ;
}

const vector<char>& PDB_t::getSSETypes(
  const uint32_t modelID,
  const uint32_t chainID
) const {
	return model[modelID].getSSETypes(chainID) ;
}

char PDB_t::getSSType(
  const uint32_t modelID, 
  const uint32_t chainID,
  const uint32_t resID
) const{
	 
	return	model[modelID].getSSType(chainID, resID);

}

void PDB_t::setSSType( 
  const uint32_t modelID, 
  const uint32_t chainID, 
  const uint32_t resID, 
  const char c 
) {
	model[modelID].setSSType(chainID,resID,c) ;
}

/*appends both the SSE boundary and type*/
const uint32_t PDB_t::getnSSEs( 
  const uint32_t modelID,
  const uint32_t chainID
) {
	return model[modelID].getnSSEs( chainID) ;
}
void PDB_t::setSSE( 
  const uint32_t modelID,
  const uint32_t chainID,
  const uint32_t s, 
  const uint32_t e, 
  const char c 
) {
	model[modelID].setSSE(chainID,s,e,c) ;
}

void PDB_t::sortSSE() {
	const uint32_t nM = getnModels() ;	
	for( uint32_t m = 0 ; m < nM ; m++ ) {
		model[m].sortSSE() ;
	}
}

void PDB_t::genSSEPymolCmds(string outdir ) {
	for( uint32_t i = 0 ; i < nModels ; i ++  ) {
		uint32_t found = fname.find_last_of("/\\") ;
		string fn = outdir+fname.substr(found+1) ;
		if( nModels == 1 ) {
			fn += ".pml" ;
		}
		else {
			fn += "_model-" ;
			char mstr[100] ;
			sprintf( mstr, "%u", i ) ;
			fn += mstr ;
			fn += ".pml" ;
		}

		ofstream out( fn.c_str(), ios::out ) ;
		assert(out ) ;

		out << "load " << fname << endl ;
		out << "hide" << endl ;
		out << "alter all, ss=\'L\'" << endl ;
		out << "show cartoon" << endl ;

		//cartoon parameters
		out << "set cartoon_fancy_helices,1" << endl ;
		out << "set cartoon_discrete_colors,1"  << endl ;
		out << "set cartoon_highlight_color, grey60" << endl ;
		out << "set cartoon_dumbbell_length,1.0" << endl ;
		out << "set cartoon_rect_length,1.10000" << endl ;
		out << "set cartoon_loop_radius,0.3" << endl ;
		out << "set cartoon_smooth_loops=0" << endl ;

		model[i].genSSEPymolCmds(out) ;

		out << "deselect" << endl ;
		out << "rebuild" << endl  ;
		out.close() ;
	}
}

/* This function generates pymol commands that color a protein 
   based on arbitrary segmentation/delineation provided to this
   routine. Before invoking this function, ensure that the segmentation
   does not have residue ranges that are overflows
*/

void PDB_t::genSegmentPymolCmds(
  vector<vector<uint32_t> > &seg, 
  vector<char> &labels, 
  vector<string> &colors,
  uint32_t modelIndx,
  uint32_t chainIndx,
  string outdir 
) {
	assert( modelIndx < nModels ) ;
	uint32_t found = fname.find_last_of("/\\") ;
	string fn = outdir+fname.substr(found+1) ;
	char mstr[100] ;
	sprintf( mstr, "_m%lu_c%lu", (unsigned long)modelIndx, (unsigned long)chainIndx ) ;
	fn += mstr ;
	fn += "_Seg.pml" ;

	ofstream out( fn.c_str(), ios::out ) ;
	assert(out ) ;

	out << "load " << fname << endl ;
	out << "hide" << endl ;
	out << "alter all, ss=\'L\'" << endl ;
	out << "show cartoon" << endl ;

	//cartoon parameters
	out << "set cartoon_fancy_helices,1" << endl ;
	out << "set cartoon_discrete_colors,1"  << endl ;
	out << "set cartoon_highlight_color, grey60" << endl ;
	out << "set cartoon_dumbbell_length,1.0" << endl ;
	out << "set cartoon_rect_length,1.40000" << endl ;
	out << "set cartoon_loop_radius,0.3" << endl ;
	out << "set cartoon_smooth_loops=0" << endl ;


	model[modelIndx].genSegmentPymolCmds(seg, labels,colors, chainIndx, out) ;

	out << "deselect" << endl ;
	out << "rebuild" << endl  ;
	out.close() ;
}

void PDB_t::genSparseRepPDB( string outdir ) {
	for( uint32_t i = 0 ; i < nModels ; i ++  ) {
		uint32_t found = fname.find_last_of("/\\") ;
		string fn = outdir+fname.substr(found+1) ;
		if( nModels == 1 ) {
			fn += "_sparse.pdb" ;
		}
		else {
			fn += "_model-" ;
			char mstr[100] ;
			sprintf( mstr, "%u", i ) ;
			fn += mstr ;
			fn += "_sparse.pdb" ;
		}

		ofstream out( fn.c_str(), ios::out ) ;
		assert(out ) ;

		model[i].genSparseRepPDB(out) ;
		out.close() ;
	}
}


/* function to find axes and center of mass*/
void PDB_t::computeAxesAndCentersOfMass() {
	const uint32_t nM = getnModels() ;	
	for( uint32_t m = 0 ; m < nM ; m++ ) {
		model[m].computeAxesAndCentersOfMass() ;
	}
}


const vector<vector<double> >& PDB_t::getSSEAxes(
  const uint32_t modelID,
  const uint32_t chainID
) const {
	return model[modelID].getSSEAxes(chainID) ;
}

const vector<vector<double> >& PDB_t::getSSECentersOfMass(
  const uint32_t modelID,
  const uint32_t chainID
) const {
	return model[modelID].getSSECentersOfMass(chainID) ;
}

const vector<vector<double> >& PDB_t::getSSEStartPoints(
  const uint32_t modelID,
  const uint32_t chainID
) const {
	return model[modelID].getSSEStartPoints(chainID) ;
}

const vector<vector<double> >& PDB_t::getSSEEndPoints(
  const uint32_t modelID,
  const uint32_t chainID
) const {
	return model[modelID].getSSEEndPoints(chainID) ;
}



/* tableau related function*/
void PDB_t::computeTableau() {
	const uint32_t nM = getnModels() ;	
	for( uint32_t m = 0 ; m < nM ; m++ ) {
		model[m].computeTableau() ;
	}
}

const vector<vector<double> >& PDB_t::getTableau_orientations(
  const uint32_t modelID,
  const uint32_t chainID
) const {
	return model[modelID].getTableau_orientations(chainID) ;
}

const vector<vector<uint8_t> >& PDB_t::getTableau_codes(
  const uint32_t modelID,
  const uint32_t chainID
) const {
	return model[modelID].getTableau_codes(chainID) ;
}

void PDB_t::printTableauInfo() {
	uint32_t found = fname.find_last_of("/\\");
  	string fn = fname.substr(found+1)+".tableau_info" ;

	ofstream out( fn.c_str(), ios::out ) ;
	assert(out ) ;
	for( uint32_t i = 0 ; i < nModels ; i ++  ) {
		out << "----- MODEL " << i << " -----" << endl ;
		model[i].printTableauInfo(out) ;
	}
	out.close() ;
}


/*
int main( int argc, char *argv[] ) {
//	PDB_t pdb( "./data/1HHO.pdb" );
//	PDB_t pdb( "./data/1a62_MSEexample.ent" );
//	PDB_t pdb( "./data/3myv_anisouexample.pdb" );
//	PDB_t pdb( "./data/1h4w_insertioncodeexample.pdb" );
//	PDB_t pdb( "./data/2kj3_multimodelexample.pdb" );
//	pdb.print() ;

	PDB_t pdb( argv[1] ) ;
	uint32_t nModels, nChains, nRes, nAtoms ;
	cout << "nModels = " << (nModels=pdb.getnModels()) << endl ;
	//pdb.print() ;

	for( uint32_t i = 0 ; i < nModels ; i ++ ) {
		cout << "nChains in Model-" << i << " = " 
			<< (nChains = pdb.getnChains( i )) << endl ;
		for( uint32_t j = 0 ;  j < nChains ; j ++ ) {
			cout << "\tnRes in Chain-" << j << " = " 
			<< (nRes = pdb.getnResidues( i, j )) << endl  ;
			for( uint32_t k = 0 ;  k < nRes ; k++ ) {
				cout << "\t\tnAtoms in Res-" << k << " = " 
				<< (nAtoms = pdb.getnAtoms( i, j, k ))  ;
				cout << "\tRescode ("
				<< pdb.getResCode( i, j,k ) <<")" << endl ;
				//cout << "\tSS type" << pdb.getSSType(i,j,k) << endl ;
				const double *coords ;
				char patt[] = "CA" ;
				ATOM_t natom = 
				   pdb.getAtomRecord(i,j,k,patt) ;
				if( natom.isEmpty() ) {
					cout << "Does not have a atom:"
					       <<  patt << endl ;	
				}
				else {
					coords = natom.getXYZ() ;
					cout << "\t\t\t" << coords[0] << " " ;
					cout <<  coords[1] << " " ;
					cout <<  coords[2] << endl ;
					natom.print() ;
				}
			}
		}
	}
	return 0 ;
}
*/
