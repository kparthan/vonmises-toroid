#ifndef PDB_H
#define PDB_H
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <string>
#include <sstream>
#include <stdint.h>
#include <map>
#include <vector>
#include <stdint.h>
#include <math.h>
#include "geometry3D.h"
#include "superpose3D.h"

using namespace std;

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
* 18-20 | 17-19 | Residue name|  resiName | Residue name.             *
* 22    | 21    | Character   |  chainID  | Chain identifier.         *
* 23-26 | 22-25 | Integer     | resiSeqNum| Residue sequence number.  *
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

/*#### ATOM_t CLASS ###*/
class ATOM_t {
	private:
		string verbatimLine ;
		string recordType ;
		uint32_t serial ;
		string atomName ;
		char altLoc ;
		string resiName ;
		char chainID ;
		int32_t resiSeqNum ;
		char iCode ;
		double coord[3] ;
		double x, y, z ;
		double occupancy ;
		double tempFactor ;
		//string element ;
		//string charge ;
		bool EMPTY ;
	public:
		ATOM_t() ;
		void process( const char [] ) ;	
		const double* getXYZ() const;
		string getResIdentifier() const ;
		 int32_t getResSeqNum() const ;
		double* getCoords() ; // same as getXYZ but returns non const
		void print() ;
		void printPDB(ofstream&) ;
		void printPDB(ofstream&, Superpose3DClass&) ;
		bool isEmpty() ;
};

/*#### RESIDUE_t CLASS ###*/
class RESIDUE_t {
	private: 
		uint32_t nAtoms ;
		char resType[4] ; // 3 letter code
		char resCode ;
		map<string,int> atomName_indx;
		vector< ATOM_t > atom ;
		double dihedralsPhiPsiOmega[3] ; // phi psi, omega
		
		/*stores the residue index that forms an Hbond with this residue
		 [0] stores the index of res which accepts  the Hbond (N-H --> O=C)
		 [1] stores the index of res which receives the Hbond (C=O <-- H=N)
		 */
		uint32_t NHCO_HBondResIndx[2] ;  

		/* Variables related to Secondary Structure calculations.
		 * The following variables can be used only after the
		 * PDB coordinates are parsed and the PDB object whose
		 * secondary structural definitions have to be computed
		 * is passed through an instantiation of an independent 
		 * SSE_t class.
		*/
		// Note the following symbols are used:
		// 'C' Cis-Peptide
		// 'G' 3_{10} helix
		// 'H' Alpha helix
		// 'I' Pi helix
		// 'E' beta strand.
		// 'L' for a sheet
		// '-' unassigned. This suggests that the Secondary
		//       stucture calculation was not yet performed.
		char ssType ; // this is a residue level ss definition.
		bool ssDefFlag ; // flag that tells you if SS is defined 
	public:
		RESIDUE_t() ;	
		const uint32_t getnAtoms() const ;
		const ATOM_t& getAtomRecord( const char[] ) ;
		double* getCoords( const char[] ) ; 
		double* getCACoords() ;
		char getResCode() const ;
		string getResIdentifier() const ;
		 int32_t getResSeqNum() const ;
		void process( const char [] ) ;
		void aaAssign( const char [], const char ) ;
		bool isAtomIndexed( const string ) ;
		void setPhiPsiOmega( const double, const double, const double) ;
		void setNHCO_HBondResIndx( const uint32_t, const uint32_t) ;
		uint32_t* getNHCO_HBondResIndx() ;
		double* getPhiPsiOmega();
		void print() ;
		void printPDB(ofstream&) ;
		void printPDB(ofstream&, Superpose3DClass&) ;

		//ss related member functions.
		void setSSType( const char ) ; 
		char getSSType() const ;
};

/*#### CHAIN_t CLASS ###*/
class CHAIN_t {
	private:
		char chainID ;
		uint32_t nResidues ;
		//While resSeq is an integer aa id, it alone is not
		//enough to disambiguate amino acids in the chain.
		//This is due to the use of insertion code (iCode)
		//to have a numbering that allows two proteins to be
		//comparable. Therefore the map below uses resSeq+iCode
		//to identify uniquely various amino acids.
		map<string,int> resNum_indx ;
		map<string,char> aacode ;
		vector<RESIDUE_t> residue ;

		/* Variables related to Secondary Structure calculations.
		 * These will be filled by methods from SSE_t class. */
		// data structure to store Start and End points of 
		// Secondary structural elements (SSEs). Here SSEs 
		// considered are Helices and strands and nothing else. 
		// No distinction is made between 3_{10}, Alpha and Pi helix.
		// For residue level definition, use the methods and vars
		// in residues objects corresponding to the chain obj.
		uint32_t nSSEs ;
		vector<vector<uint32_t> > sse ;
		vector<char> sseType ; // an associative list to the above.
		bool sseFlag ;

		/* once the SSEs are computed, an axis for each SSE
		 * and CM is stored in the data structures below.*/
		vector<vector<double> > axes ;  
		vector<vector<double> > centersOfMass ;
		/* In addition, the first and last CA coordinates of each
		 * SSE is projected on to the SSE's axis and start and end
		 * points are computed. Note that center of Mass computed 
		 * above need NOT be the same as the mid point between
		 * start and end points computed here. */
		vector<vector<double> > startPoints ;
		vector<vector<double> > endPoints ;

		bool axesAndCMFlag ;

		/*once the Axes and CM are computed, tableau related
		 *data is stored in the data structures below*/
		bool tableauFlag ;
		//orientation angles of all pairs of SSEs
		vector<vector<double> > omega ; 
		//tableau code:
		//PE=0 //PD=1 //RD=2 //RT=3
		//OT=4 //OS=5 //LS=6 //LE=7
		vector<vector<uint8_t> > tableau ;
		// converts an orientation angle to a uint8_t code.
		// An argument uint8_t has to be passed which
		// will contain the on return
		void omegaToBytecode( const double, uint8_t& ) ;

		// converts an orientation angle to a double quadrant
		// encoding.
		// An argument char[2] has to be passed which will be 
		// filled in with the double character encoding of a 
		// tableau
		void omegaToEncoding( const double, char[2] ) ;

		void byteCodeToEncoding( const uint8_t, char[2] ) ;
		
		//colors related var for generating pymol load script
		uint32_t nPymolColors ;
		vector<string> pymolColor ;
		void setPymolColors() ;
		bool computeAxisAndCM( const uint32_t ) ;

	public:
		CHAIN_t() ;
		void process( const char [] ) ;	
		void setChainID( const char ) ;
		const uint32_t getnResidues() const ;
		const uint32_t getnAtoms( const uint32_t ) const ; 
		char getResCode( const uint32_t ) const ; 
		string getResIdentifier( const uint32_t ) const ; 
		 int32_t getResSeqNum( const uint32_t) const ;
		const ATOM_t& getAtomRecord( const uint32_t, const char[] ) ; 
		double* getCACoords( const uint32_t ) ; 
		char getChainID() ;
		void computePhiPsiOmega() ;
		void compute_NHCO_HBond() ;
		double* getPhiPsiOmega( const uint32_t);
		uint32_t* getNHCO_HBondResIndx( const uint32_t);
		bool isResidueIndexed( const char[] ) ;
		void print() ;
		void printPDB(ofstream&) ;
		void printPDB(ofstream&, Superpose3DClass&) ;


    string getAminoAcidSequenceOfChain();
    vector<string> getResIdentifiersOfChain();
    vector<vector<double> > getCACoordsOfChain();
    void getCACoordsOfChain(vector<vector<double> >&, string&, vector<string>&, vector<uint32_t>&);
    
		//sse related member functions
		const uint32_t getnSSEs() const  ;
		const vector<vector<uint32_t> >& getSSEBoundaries() const ; 
		const vector<char>& getSSETypes()const  ;
		char getSSType( const uint32_t ) const ;

		void setSSType( const uint32_t, const char ) ;
		void setSSE(const uint32_t , const uint32_t, const char ) ;
		//sorts SSE boundaries and associated types.
		void sortSSE() ;
		void genSSEPymolCmds(ofstream &) ; 
		void genSegmentPymolCmds(vector<vector<uint32_t> >&,vector<char>&, vector<string>&, ofstream &) ; 
		void genSparseRepPDB(ofstream &) ; 
		void computeAxesAndCentersOfMass() ;
		const vector<vector<double> >& getSSEAxes() const ; 
		const vector<vector<double> >& getSSECentersOfMass() const; 
		const vector<vector<double> >& getSSEStartPoints() const; 
		const vector<vector<double> >& getSSEEndPoints() const; 
		//tableau related
		void computeTableau() ;
		const vector<vector<double> >& getTableau_orientations() const ; 
		const vector<vector<uint8_t> >& getTableau_codes() const; 
		void printTableauInfo(ofstream &) ; 
};


/*#### MODEL_t CLASS ###*/
class MODEL_t {
	private:
		uint32_t nChains ;
		uint32_t modelID ; 
		map<char,int> chainID_indx;
		vector<CHAIN_t> chain ;

	public:
		MODEL_t() ;
		void process( const char [] ) ;
		bool isChainIndexed(const char, uint32_t&) ;
		const uint32_t getnChains() const ;
		const uint32_t getnResidues( const uint32_t ) const ;
		const uint32_t getnAtoms( const uint32_t, const uint32_t ) const ;
		char getResCode( const uint32_t, const uint32_t ) const ;
		string getResIdentifier( const uint32_t, const uint32_t ) const ;
		 int32_t getResSeqNum( const uint32_t, const uint32_t ) const ;
		const ATOM_t& getAtomRecord( const uint32_t, const uint32_t, const char[] ) ;
		double* getCACoords( const uint32_t, const uint32_t ) ;
		char getChainID( const uint32_t ) ;
		void computePhiPsiOmega() ;
		void compute_NHCO_HBond() ;
		double* getPhiPsiOmega( const uint32_t, const uint32_t);
		uint32_t* getNHCO_HBondResIndx( const uint32_t, const uint32_t);
		void print() ;
		void printPDB(ofstream&) ;
		void printPDB(ofstream&, Superpose3DClass&) ;

    string getAminoAcidSequenceOfChain(const uint32_t);
    vector<string> getResIdentifiersOfChain(const uint32_t);
    vector<vector<double> > getCACoordsOfChain(const uint32_t);
    void getCACoordsOfChain(const uint32_t, vector<vector<double> >&, string&, vector<string>&, vector<uint32_t>&);
    void getCACoordsOfChain(const char, vector<vector<double> >&, string&, vector<string>&, vector<uint32_t>&);

		//secondary structure related functions
		const uint32_t getnSSEs( const uint32_t ) ;
		const vector<vector<uint32_t> >& getSSEBoundaries( const uint32_t ) const ;
		const vector<char>& getSSETypes( const uint32_t ) const ;
		char getSSType( const uint32_t, const uint32_t ) const ;
		void setSSType( const uint32_t,const uint32_t,const char);
		void setSSE(const uint32_t,const uint32_t,const uint32_t,char) ;
		//sorts SSE boundaries and associated types.
		void sortSSE() ;
		void genSSEPymolCmds( ofstream &) ; 
		void genSegmentPymolCmds(vector<vector<uint32_t> >&,vector<char>&, vector<string>&, uint32_t, ofstream &) ; 
		void genSparseRepPDB( ofstream &) ; 
		void computeAxesAndCentersOfMass() ;
		const vector<vector<double> >& getSSEAxes(const uint32_t) const ; 
		const vector<vector<double> >& getSSECentersOfMass( const uint32_t) const; 
		const vector<vector<double> >& getSSEStartPoints( const uint32_t) const; 
		const vector<vector<double> >& getSSEEndPoints( const uint32_t) const; 

		//tableau related
		void computeTableau() ;
		const vector<vector<double> >& getTableau_orientations(const uint32_t) const ; 
		const vector<vector<uint8_t> >& getTableau_codes(const uint32_t) const; 
		void printTableauInfo(ofstream &) ; 
};


/*### PDB_t CLASS ###*/
class PDB_t {
	private: 
		string fname ;
		string title ;// stores title record if it exists
		uint32_t nModels ;
		vector< MODEL_t > model ;  
		bool parsePDB( const string ) ;
		bool isValidRecord( const char[] ) ;
	public:
		PDB_t(const string fname) ;
		// Member functions to access summary/properties of the PDB file.
      const string getTitle() const ;
		const string getTitleSummary() const ;
		string getFileName() const ;
		const uint32_t getnModels() const ;
		const uint32_t getnChains( const uint32_t ) const;
		const uint32_t getnResidues( const uint32_t, const uint32_t ) const ;
		const uint32_t getnAtoms( const uint32_t, const uint32_t, const uint32_t ) const ;
      // Member functions to access information at a SPECIFIED RESIDUE level.
      // They all will need THREE integers: modelIndx, chainIndx, and
      // residueIndx to access residue level infomation.
		char getResCode( const uint32_t, const uint32_t, const uint32_t ) const ;
		string getResIdentifier( const uint32_t, const uint32_t, const uint32_t ) const ;
		 int32_t getResSeqNum( const uint32_t, const uint32_t, const uint32_t ) const ;
		const ATOM_t& getAtomRecord( const uint32_t, const uint32_t, const uint32_t, const char []);
		double* getCACoords( const uint32_t, const uint32_t, const uint32_t );
      string getAllChainIDs(const uint32_t);

      // Member functions to access information at a SPECIFIED CHAIN level.
      // They all will need TWO integers: modelIndx and chainIndx to access 
      // chain level information.
		char getChainID( const uint32_t, const uint32_t) ;
      string getAminoAcidSequenceOfChain(const uint32_t, const uint32_t);
      vector<string> getResIdentifiersOfChain(const uint32_t, const uint32_t);
      vector<vector<double> > getCACoordsOfChain(const uint32_t, const uint32_t);
      void getCACoordsOfChain(const uint32_t, const uint32_t, vector<vector<double> >&, string&, vector<string>&, vector<uint32_t>&);
      void getCACoordsOfChain(const uint32_t, const char, vector<vector<double> >&, string&, vector<string>&, vector<uint32_t>&);
      void getCACoordsOfChains(const uint32_t, const string, vector<vector<vector<double> > >&, vector<string> &, vector<vector<string> >&, vector<vector<uint32_t> >&);
    
		void computePhiPsiOmega() ;
		void compute_NHCO_HBond() ;
		double* getPhiPsiOmega(const uint32_t, const uint32_t, const uint32_t);
		uint32_t * getNHCO_HBondResIndx(const uint32_t, const uint32_t, const uint32_t);

		void print() ;
		void printPDB(string) ;
		void printPDB(string, Superpose3DClass&) ;

		//secondary structure related functions
		const uint32_t getnSSEs( const uint32_t , const uint32_t) ;
		const vector<vector<uint32_t> >& getSSEBoundaries( 
		const uint32_t, const uint32_t ) const ;
		const vector<char>& getSSETypes( 
		const uint32_t, const uint32_t ) const ;
		char getSSType( 
		const uint32_t, const uint32_t, const uint32_t ) const ;

		void setSSType(const uint32_t,const uint32_t,const uint32_t,const char);
		void setSSE(const uint32_t,const uint32_t,const uint32_t,const uint32_t,char) ;
		//sorts SSE boundaries and associated types.
		void sortSSE() ;
		void genSSEPymolCmds( string ) ;
		void genSegmentPymolCmds(vector<vector<uint32_t> >&, vector<char>&, vector<string>&, uint32_t, uint32_t, string ) ; 
		void genSparseRepPDB( string ) ;
		void computeAxesAndCentersOfMass() ;
		const vector<vector<double> >& getSSEAxes(const uint32_t, const uint32_t) const ; 
		const vector<vector<double> >& getSSECentersOfMass( const uint32_t, const uint32_t) const; 
		const vector<vector<double> >& getSSEStartPoints( const uint32_t, const uint32_t) const; 
		const vector<vector<double> >& getSSEEndPoints( const uint32_t, const uint32_t) const; 

		//tableau related
		void computeTableau() ;
		const vector<vector<double> >& getTableau_orientations(const uint32_t, const uint32_t) const ; 
		const vector<vector<uint8_t> >& getTableau_codes( const uint32_t, const uint32_t) const; 
		void printTableauInfo() ; 
} ;

extern ATOM_t tmp1 ;

#endif
