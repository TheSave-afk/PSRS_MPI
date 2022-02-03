using namespace std;
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <string.h>
#include <vector>
#include <queue>
#include <math.h>
#include "mpi.h"

int compara(const void *i, const void *j)
{
	int int1 = *reinterpret_cast<const int *>(i);
	int int2 = *reinterpret_cast<const int *>(j);
	if (int1<int2) return -1;
	if (int1>int2) return 1;
	return 0;
}

// funzione che verifica se l'array è ordinato
string checkSort(int array[], int arrayStart, int arrayLength)
{
	for(int i=arrayStart; i<arrayStart+arrayLength-1;i++)
	{
		if (array[i]>array[i+1])
			return "non è ";
	}

	return " è ";
}

string checkSort(int array[], int arrayLength)
{
	return checkSort(array,0,arrayLength);
}

struct m_data
{
	//index dell'array da unire
	int st_index;
	//index all'elemento nell'array
	int index;
	//valore dell'elemento
	int st_value;

	m_data(int st=0, int id=0, int stv = 0):st_index(st),index(id),st_value(stv){}

};


// comparison operator
bool operator<( const m_data & uno, const m_data & due)
{
	return uno.st_value > due.st_value;
}

/*
 * Un array bidimensionale di regularSamples per ogni processo, la lunghezza dei regularSamples
 * per ogni processo, il numero di array da unire, l'array dei risultati, il volume totale di dati da unire
 */ 
int multiple_merge(int * starts[], const int len[], const int n,int newArray[], const int newArray_Length)
{

 	priority_queue< m_data> priorities;
	/* 
     *   Aggiungi il primo numero di ciascun array alla coda di priorità 
	 * 	 se la lista non è vuota
     */
 	for(int i = 0; i < n ; i++)
 	{
		if (len[i] > 0)
		{
			priorities.push(m_data(i,0,starts[i][0]));
		}
	}

	int newArray_index = 0; 
	while (!priorities.empty() && (newArray_index<newArray_Length))
	{
		// seleziono l'elemento più piccolo e lo toglie dalla coda prioritaria
		m_data array = priorities.top();
		priorities.pop();

		// inserisco l'elemento più piccolo nell'array ordinato
		newArray[newArray_index++] = starts[array.st_index][array.index];

		// se i dati acquisiti non sono l'ulitmo elemento nell'array allora
		// inserisco nella coda prioritaria l'elemento successivo
		if ( len[array.st_index] > (array.index + 1))
		{
			priorities.push(m_data(array.st_index, array.index + 1,
								starts[array.st_index][array.index + 1]));
		}
}

return newArray_index;
}

int main( int argc, char* argv[])
{
	// ID del processore, numero totale di processori
    int ID, nProc;

    double startwtime = 0.0, endwtime;

	// inizializzazione dell'ambiente di esecuzione MPI
    MPI::Init(argc, argv);

	//numero di processori e ID
    nProc = MPI::COMM_WORLD.Get_size();
    ID = MPI::COMM_WORLD.Get_rank();

	// -DS usa myData_Size, ovvero il numero di campioni passati da riga di comando
	// -SR usa seed_rand
	int seed_rand = 1000;
	int myData_Size = 4000;
	for(int i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-DS")==0)
		{
			myData_Size = atoi(argv[i+1]); i++;
		}
		else if (strcmp(argv[i],"-SR")==0)
		{
			seed_rand = atoi(argv[i+1]); i++;
		}
	}

    // array per immagazzinare i dati iniziali non ordinati
	// myData_Lengths[] e myData_Starts[] usati con Scatterv() per distribuire
	// myData[] agli altri processori
	int myData[myData_Size];
	int myData_Lengths[nProc];
	int myData_Starts[nProc];

	// buffer di comunicazione usato per la determinazione dei valori dei pivot
	int pivot_buff[nProc*nProc];
	int pivot_buffSize;


	// Calcolo la lunghezza di mydata[] per distribuirla
	// al processore nProc. L'ultimo processore avrà lunghezza extra
	for(int i=0;i<nProc;i++)
	{
		myData_Lengths[i] = myData_Size/nProc;
		myData_Starts[i]= i*myData_Size/nProc;
	}
	myData_Lengths[nProc-1]+=(myData_Size%nProc);


	// Il nodo principale inizializza i campioni e fa partire il tempo di esecuzione
    if (ID == 0)
    {
		// seme random
		srandom(seed_rand);
		for(int index=0; index<myData_Size; index++)
		{
			myData[index] = random()% 900;
		}

		startwtime = MPI::Wtime();
    }

	/*
	 * il processorore principale distribuisce i dati a tutti i processori
	 * ed esegue in locale
	 */
	if (ID == 0)
	{
		MPI::COMM_WORLD.Scatterv(myData,myData_Lengths,myData_Starts,MPI::INT,MPI_IN_PLACE,myData_Lengths[ID],MPI::INT,0);
	}
	else
	{
		MPI::COMM_WORLD.Scatterv(myData,myData_Lengths,myData_Starts,MPI::INT,myData,myData_Lengths[ID],MPI::INT,0);
	}

	//Ogni processo ordina i propri dati
	qsort(myData,myData_Lengths[ID], sizeof(int), compara);


	//Tutti i processi raccolgono campioni regolari dalla lista ordinata con offset = index*myData
	for(int index=0;index<nProc;index++)
	{
		pivot_buff[index]= myData[index*myData_Lengths[ID]/nProc];
	}

	//Il processore con rank 0 raccoglie i pivot, i suoi sono in locale
	if (ID == 0)
	{
		MPI::COMM_WORLD.Gather(MPI_IN_PLACE,nProc,MPI::INT,
			pivot_buff,nProc,MPI::INT,0);
	}
	else
	{
		MPI::COMM_WORLD.Gather(pivot_buff,nProc,MPI::INT,
			pivot_buff,nProc,MPI::INT,0);
	}

	if (ID == 0)
	{
		
		int *starts[nProc];  //array delle liste
		int lengths[nProc];  //array della lunghezza delle liste

		//unisce le liste in una sola
		for(int i = 0; i < nProc; i++)
		{
			starts[i]=&pivot_buff[i*nProc];
			lengths[i]=nProc;
		}
		int tempbuffer[nProc*nProc];
		multiple_merge(starts,lengths,nProc,tempbuffer,nProc*nProc);

		//seleziona nProc-1 pivot da inviare
		for(int i=0; i<nProc-1; i++)
		{
			pivot_buff[i] = tempbuffer[(i+1)*nProc];
		}
	}
	// broadcast dei valori della partizione
	MPI::COMM_WORLD.Bcast(pivot_buff,nProc-1,MPI::INT,0);


	int class_Start[nProc];  //indice dell'inizio della i-ma classe
	int class_Length[nProc]; //lunghezza della i-ma classe 


	// ogni processo partiziona la lista usando il valore di pivot_buff[]
	int dataindex=0;
	for(int classindex=0; classindex<nProc-1; classindex++)
	{
		class_Start[classindex] = dataindex;
		class_Length[classindex]=0;

		// dataindex si riferisce ai dati nella classe corrente
		while((dataindex< myData_Lengths[ID])
			&& (myData[dataindex]<=pivot_buff[classindex]))
		{
			class_Length[classindex]++;
			dataindex++;
		}
	}
	// imposta start and lenght per l'ultima classe
	class_Start[nProc-1] = dataindex;
	class_Length[nProc-1] = myData_Lengths[ID] - dataindex;

	//Tutte le i-esime classi sono raccolte dal processore i
	int recv_buff[myData_Size];    //buffer per tenere in memoria tutti i membri della classe i
	int recv_Lengths[nProc];     // numero del processore, lunghezza di ogni numerodiprocessore-iesima classe
	int recv_Starts[nProc];      //indice che ci dice dove iniziare a memorizzare


	//si prendono tutte le liste orindate e tutti i valori ordinati del iprocessore-esima classe
	for(int iprocessor=0; iprocessor < nProc; iprocessor++)
	{
		//ogni processo raccoglie la lunghezza degli valori ordinati nella classe dell' iprocessor
		MPI::COMM_WORLD.Gather(&class_Length[iprocessor], 1, MPI::INT,recv_Lengths,1,MPI::INT,iprocessor);

		// da questa lunghezza sono calcolati gli offset di start su ogni processore
		if (ID == iprocessor)
		{
			recv_Starts[0]=0;
			for(int i=1;i<nProc; i++)
			{
				recv_Starts[i] = recv_Starts[i-1]+recv_Lengths[i-1];
			}
		}

		MPI::COMM_WORLD.Gatherv(&myData[class_Start[iprocessor]],class_Length[iprocessor],MPI::INT,
			recv_buff,recv_Lengths,recv_Starts,MPI::INT,iprocessor);
	}


	// multiple_merge
	int *m_Starts[nProc];
	for(int i=0;i<nProc;i++)
	{
		m_Starts[i]=recv_buff+recv_Starts[i];
	}
	multiple_merge(m_Starts,recv_Lengths,nProc,myData,myData_Size);

	int mysend_Length = recv_Starts[nProc-1] + recv_Lengths[nProc-1];

	int send_Lengths[nProc]; // lunghezza delle classi
	int send_Starts[nProc];  // starting points delle classi
	// il processo con rank 0 raccoglie le lunghezze di tutti i dati
	MPI::COMM_WORLD.Gather(&mysend_Length,1,MPI::INT,send_Lengths,1,MPI::INT,0);

	// il nodo con rank 0 calcola l'offset di star dalla lunghezza delle classi da raccogliere
	if (ID == 0)
	{
		send_Starts[0]=0;
		for(int i=1; i<nProc; i++)
		{
			send_Starts[i] = send_Starts[i-1]+send_Lengths[i-1];
		}
	}

	// il processo con rank 0 raccoglie le liste e le ordina
	int arraySorted[myData_Size];
	MPI::COMM_WORLD.Gatherv(myData,mysend_Length,MPI::INT,arraySorted,send_Lengths,send_Starts,MPI::INT,0);


    double tseq = myData_Size * double(log10(myData_Size)) ;
    double speed_up = tseq / -(endwtime-startwtime) ;
    double efficienza = speed_up / myData_Size;

    if (ID == 0)
	{
		endwtime = MPI::Wtime();
        cout << "\nTempo di esecuzione totale (in secondi) = "
		     <<  setprecision(4) << endwtime-startwtime << endl;

		cout << "\nESITO TEST: il data set " << checkSort(arraySorted,myData_Size) << " ordinato!" << endl;

		cout << "\nPRESTAZIONI:\n"
		     <<  "\tSPEED-UP = " << speed_up << "\n\tEFFICIENZA = " << efficienza <<"\n"<< endl;
	}

    // chiusura dell'ambiente MPI
    MPI::Finalize();
    return 0;
}
