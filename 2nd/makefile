SOURCES=c.c \
	vector.c\
	runge.c

OBJDIR=OBJECTS

PROG_NAME=a.out

CFLAGS=-O3 -g -lm -Wall -Wextra -Werror

all: $(PROG_NAME)


OBJS=$(SOURCES:%.c=$(OBJDIR)/%.o)




clean:
	rm -rf $(OBJDIR) $(PROG_NAME) data.res
Clean:
	rm *~


ifneq ($(MAKECMDGOALS),clean)
-include $(SOURCES:%.c=$(OBJDIR)/%.d)
endif

$(OBJDIR)/%.d: %.c
	mkdir -p $(OBJDIR) $(OBJDIR)/rooms
	$(CC) -MM -MF $@ -E $< -MT '$@ $(OBJDIR)/$*.o' $(CFLAGS)

$(OBJDIR)/%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(PROG_NAME): $(OBJS)
	$(CC)  -s -o $@ $+ $(CFLAGS)
