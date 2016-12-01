CC=gcc
SRC=src/
EXE1=potencia
EXE2=primo

todo: $(EXE1) $(EXE2)

$(EXE1): $(SRC)$(EXE1).c
	$(CC) -g -o $@ $< -lgmp -lm
	
$(EXE2): $(SRC)$(EXE2).c
	$(CC) -g -o $@ $< -lgmp -lm
