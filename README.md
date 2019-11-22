# MscThesis

In this repository you can find the codes of the programs developed for my thesis work done at Laboratoire de Physique Théorique et Modèles statistiques of Université Paris-Sud (France), during my MSc in Physics at Università degli Studi di Padova (Italy).

- The program "pavoidfree2D.py" generates a desired number of free directed polymers in 1+1 dimensions and counts which avoid each other and which cross. Setting the size of the polymers `T`, their initial (and final) distance `dx` and the number of samples `trials`, the program produces as output an estimation of the non-crossing probability, called `pavoid` for that size and initial distance.

- The program "pavoidfree3D.cpp" generates a desired number of free directed polymers in 2+1 dimensions and counts which avoid each other and which cross. Setting the size of the polymers `T`, their initial (and final) distance `dx` and the number of samples `trials`, the program produces as output a .txt file of values whose arithmetical average gives an estimation of the non-crossing probability for that size and initial distance.

- The program "pavoidnoise2D.cpp" generates a desired number of ground state directed polymers in disordered media in 1+1 dimensions, thanks to a variation of the Dijkstra's algorithm, and counts which avoid each other and which cross. Setting the size of the polymers `T`, their initial (and final) distance `dx` and the number of samples `trials`, the program produces as output a .txt file of values whose arithmetical average gives an estimation of the non-crossing probability for that size and initial distance.

- The program "pavoidnoise3D.cpp" generates a desired number of ground state directed polymers in disordered media in 2+1 dimensions, thanks to a variation of the Dijkstra's algorithm, and counts which avoid each other and which cross. Setting the size of the polymers `T`, their initial (and final) distance `dx` and the number of samples `trials`, the program produces as output a .txt file of values whose arithmetical average gives an estimation of the non-crossing probability for that size and initial distance.

- The program "energygap2D.cpp" generates a desired number of ground state + Ist excited in 2 dimensions thanks to a variation of the Dijkstra's algorithm, and calculate their energy gaps. Setting the size of the polymers `T` and the number of samples `trials`, the program produces as output a .txt file containing the energy gaps exstimated for every couple of polymers.

- The program "energygap3D.cpp" generates a desired number of ground state + Ist excited in 2 dimensions thanks to a variation of the Dijkstra's algorithm, and calculate their energy gaps. Setting the size of the polymers `T` and the number of samples `trials`, the program produces as output a .txt file containing the energy gaps exstimated for every couple of polymers.

- The program "crosstimefree2D.cpp" generates a desired number of free directed polymers in 2+1 dimensions and counts the time at which the first cross occurs. Setting the size of the polymers `T`, their initial (and final) distance `dx` and the number of samples `trials`, the program produces as output a .txt file of values of the first cross time for every couple of polymer.

## How to make every C++ code work

- Compile the code 'codename.cpp':
g++ -std=c++14 codename.cpp exename
- Make the cod run:
./exename

## How to make every Python code work

- Make the code 'codename.py' run:
python codename.py

For any questions mail me at federico.lanza.5@gmail.com
