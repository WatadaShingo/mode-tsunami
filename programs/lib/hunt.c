void hunt(xx,n,x,jlo)
int n,*jlo;
float xx[],x;
{
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n-1] > xx[0]);
	if (*jlo <= -1 || *jlo > n-1) {
		*jlo=-1;
		jhi=n;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n-1) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n-1) {
					jhi=n;
					break;
				}
			}
		} else {
			if (*jlo == 0) {
				*jlo=-1;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 0) {
					*jlo=-1;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}
void huntd(xx,n,x,jlo)
int n,*jlo;
double xx[],x;
{
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n-1] > xx[0]);
	if (*jlo <= -1 || *jlo > n-1) {
		*jlo=-1;
		jhi=n;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n-1) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n-1) {
					jhi=n;
					break;
				}
			}
		} else {
			if (*jlo == 0) {
				*jlo=-1;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 0) {
					*jlo=-1;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}
