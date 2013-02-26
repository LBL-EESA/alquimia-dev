/* -*-  mode: c; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* header file to replace stdbool.h for windows visual studio
 * compiler */

#ifndef X_STDBOOL_H_
#define X_STDBOOL_H_

#ifndef __cplusplus

#define	bool	_Bool
typedef	int   _Bool;

#define	false 0
#define	true  1

#endif /* __cplusplus */

#endif     /* X_STDBOOL_H_ */
