OBJS 	= main.o hashtable.o list.o voter.o mvote.o prompt.o
SOURCE	= main.c hashtable.c list.c voter.c mvote.c prompt.c
HEADER	= hashtable.h voter.h list.h mvote.h prompt.h
OUT		= mvote
CC		= gcc
FLAGS	= -g -c

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) 

main.o: main.c
	$(CC) $(FLAGS) main.c

mvote.o: mvote.c
	$(CC) $(FLAGS) mvote.c

prompt.o: prompt.c
	$(CC) $(FLAGS) prompt.c

hashtable.o: hashtable.c
	$(CC) $(FLAGS) hashtable.c

list.o: list.c
	$(CC) $(FLAGS) list.c

voter.o: voter.c
	$(CC) $(FLAGS) voter.c

clean:
	rm -f $(OBJS) $(OUT)

count:
	wc $(SOURCE) $(HEADER)
