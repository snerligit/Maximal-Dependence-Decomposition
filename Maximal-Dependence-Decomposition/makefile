.phony:clean compile link debug run leakcheck help all

#TRAIN = "/Users/snerli/Desktop/Personal/SJSU/CS298/Extract_Variable_Length_Sites/EI_true.seq_9.out"
#TEST = "/Users/snerli/Desktop/Personal/SJSU/CS298/Extract_Variable_Length_Sites/EI_true.seq_9.out"

#TRAIN = "/Users/snerli/Desktop/Personal/SJSU/CS298/Datasets/EI/9-mer/EI_positives_train_1.txt"
#TEST = "/Users/snerli/Desktop/Personal/SJSU/CS298/Datasets/EI/9-mer/EI_positives_test_1.txt"
#OUTPUT = "/Users/snerli/Desktop/Personal/SJSU/CS298/Results/MDD/EI/9-mer/ss_1.txt"

#$(foreach number, 1 2 3 4 5 6 7 8 9 10, ./a.out /Users/snerli/Desktop/Personal/SJSU/CS298/Results/MDD/EI/7-mer/nss_$(number).txt $(CLASS) $(WINDOW_SIZE) $(GROUPS) /Users/snerli/Desktop/Personal/SJSU/CS298/Datasets/EI/7-mer/EI_positives_train_$(number).txt /Users/snerli/Desktop/Personal/SJSU/CS298/Datasets/EI/7-mer/EI_negatives_test_$(number).txt $(NO_TRAIN_SEQ) $(LEN_TRAIN_SEQ);)

#$(foreach number, 1 2 3 4 5 6 7 8 9 10, ./a.out /Users/snerli/Desktop/Personal/SJSU/CS298/Results/PWM/EI/7-mer/nss_$(number).txt $(CLASS) $(WINDOW_SIZE) $(GROUPS) /Users/snerli/Desktop/Personal/SJSU/CS298/Datasets/EI/7-mer/EI_positives_train_$(number).txt /Users/snerli/Desktop/Personal/SJSU/CS298/Datasets/EI/7-mer/EI_negatives_test_$(number).txt $(NO_TRAIN_SEQ) $(LEN_TRAIN_SEQ);)

# Running 10-fold validations. Type the following commands under run:

# $(foreach number, 1 2 3 4 5 6 7 8 9 10, ./a.exe C:/Santrupti/Personal/SJSU/CS298/Results/MDD/Train_Authentic_Test_BRCA1/EI/nss_$(number).txt $(CLASS) $(WINDOW_SIZE) $(GROUPS) C:/Santrupti/Personal/SJSU/CS298/Datasets/EI/9-mer/EI_negatives_train_$(number).txt C:/Santrupti/Personal/SJSU/CS298/Datasets/BRCA1/BRCA1_5_negatives.txt $(NO_TRAIN_SEQ) $(LEN_TRAIN_SEQ);)

# ./a.exe C:/Santrupti/Personal/SJSU/CS298/Results/MDD/Train_Authentic_Test_HBB/sliding_window_IE_ss.txt $(CLASS) $(WINDOW_SIZE) $(GROUPS) C:/Santrupti/Personal/SJSU/CS298/Datasets/IE_true_14.txt C:/Santrupti/Personal/SJSU/CS298/Datasets/HBB/beta_globin_14_mers.txt $(NO_TRAIN_SEQ) $(LEN_TRAIN_SEQ);

CLASS = 1

NO_TRAIN_SEQ = 3000
LEN_TRAIN_SEQ = 40
GROUPS = "ACGT"
WINDOW_SIZE = 0


clean:
	if [ -a *_pwm ]; \
	then \
		rm *_pwm; \
	fi;
	rm *.o *.out output
compile:
	gcc -c main.c free.c tree.c computation.c
link:
	gcc main.o free.o tree.o computation.o -lm
run:
	./a.out /Users/snerli/Work/SJSU/MDD/BRCA1_Train_Authentic_Test_Cryptic.txt $(CLASS) $(WINDOW_SIZE) $(GROUPS) /Users/snerli/Work/SJSU/Santrupti/Personal/SJSU/CS298/Datasets/EI_true_9.txt /Users/snerli/Work/SJSU/Santrupti/Personal/SJSU/CS298/Datasets/BRCA1/BRCA1_5_Positives.txt $(NO_TRAIN_SEQ) $(LEN_TRAIN_SEQ);
debug:
	gdb a.out
leakcheck:
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-origins=yes ./a.out
help:
	./a.out -h
all:
	make compile link run
