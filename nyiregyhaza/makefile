INC = -I/media/nagyd/SSD/Documents/Egyetem/HDS/30_Projekt/PSToolboxMy -I/media/nagyd/SSD/Documents/Egyetem/HDS/30_Projekt/eigen-3.2.10
CFLAGS = -O3

dummy: dummy_teszt.cpp
	$(RM) *.o *.a dummy
	g++ dummy_teszt.cpp $(INC) $(CFLAGS) -o dummy


nyiregyhaza: nyiregyhaza.cpp
	$(RM) *.o *.a nyiregyhaza
	g++ nyiregyhaza.cpp $(INC) $(CFLAGS) -o nyiregyhaza

clean:
	$(RM) *.o *.a dummy dummy_csaba nyiregyhaza
