
// Access.
// (C) of C  : RO
// (C) of NC : RO
//        C  : RO
//        NC : RW

// Plan 1. Follow constness of the wrapper
// (C) of C  -> C  :  LINK
// (C) of NC -> C  :  LINK
// (C) of C  -> NC :  COPY
// (C) of NC -> NC :  COPY
//        C  -> C  :  LINK
//        NC -> C  :  LINK
//        C  -> NC :  COPY
//        NC -> NC :  LINK

// Plan 2. Ignore constness of the wrapper
// (C) of C  -> C  :  LINK   (SAME)
// (C) of NC -> C  :  LINK   (SAME)
// (C) of C  -> NC :  COPY   (SAME)
// (C) of NC -> NC :  LINK   (DIFF)
//        C  -> C  :  LINK   (SAME)
//        NC -> C  :  LINK   (SAME)
//        C  -> NC :  COPY   (SAME)
//        NC -> NC :  LINK   (SAME)

// Plan 3. Not allowing const of NC
// (C) of C  -> C  :  LINK
// (C) of NC -> C  :  LINK
// (C) of C  -> NC :  NOT ALLOWED
// (C) of NC -> NC :  NOT ALLOWED
//        C  -> C  :  LINK (NOT NECCESARY)
//        NC -> C  :  LINK (NOT NECCESARY)
//        C  -> NC :  NOT ALLOWED
//        NC -> NC :  LINK
