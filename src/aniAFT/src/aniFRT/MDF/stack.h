#ifndef _STACK_H_
#define _STACK_H_

template <class Type1, class Type2> class PAIR
{
public:
    PAIR(){}
    
    PAIR(Type1 obj1, Type2 obj2)
    {
        elem1 = obj1;
        elem2 = obj2;
    }
    
    bool operator==(const PAIR& pair)
    {
        if(((elem1==pair.elem1)&&(elem2==pair.elem2)) || ((elem1==pair.elem2)&&(elem2==pair.elem1)))
            return true;
        else
            return false;
    }

    Type1 elem1;
    Type2 elem2;
};

///simple stack
template <class Type> class CStack
{
    protected:
        Type *top;
        int max_sz;
	private:
		Type* mem;
		int size;
        
	public:
        CStack()
        {
            mem=0;
            top=0;
			size=0;
            max_sz=0;
        }
		CStack(int maxsize)
		{
            mem = new Type[maxsize];
			top=mem+maxsize;
            max_sz=maxsize;
            size=0;
		}
		~CStack()
		{
            if(mem)
                delete [] mem;
		}
		bool Push(const Type &value)
		{
			if(size<max_sz)
			{
				*(--top)=value;
				++size;
				return true;
			}
			return false;
		}
        bool GetTop(Type *dst) const
		{
			if(size>0)
			{
				*dst=*top;
				return true;
			}
			return false;
		}
		bool Pop(Type *dst)
		{
			if(size>0)
			{
				*dst=*top;
				top++;
				--size;
				return true;
			}
			return false;
		}
		int Size() const
		{
			return size;
		}
		bool Empty() const
		{
			return size==0;
		}
		bool Del()
		{
			if(size>0)
			{
				top++;
				--size;
				return true;
			}
			return false;
		}
		void Clear()
		{
			top+=size;
			size=0;
		}
};

///static stack (statically allocate NumElements)
template <class Type, int NumElements> class CStaticStack : public CStack<Type>
{
	protected:
	        Type *top;
	        int max_sz;
		Type massiv[NumElements];
	public:
        CStaticStack()
		{
            top = massiv+NumElements;
            max_sz = NumElements;
		}
};

/**метод выполняет быструю сортировку массива с элементами типа T, располагая его элементы по возрастанию;
\param base - массив начальных данных, \param nmemb - размер массива, \param stack - указатель на стек для рекурсии.*/ 
template <class T> void q_sort(T *base, const size_t nmemb)
{
	CStack< PAIR<int, int> > work_stack(nmemb); //стек для организации "искусственной" рекурсии
    
    PAIR<int, int> temp;
    T x;
    int r,l,i,j;
    
    if (nmemb <= 1)
        return;
    temp.elem1 = 0;
    temp.elem2 = nmemb-1;
    work_stack.Push(temp);
    do
    {
        work_stack.Pop(&temp);
        l = temp.elem1;
        r = temp.elem2;
        do
        {
            i = l;
            j = r;
            x = base[(l+r)/2];
            do
            {
                while(base[i]<x)
                    i++;
                while(x<base[j])
                    j--;
                if(i<=j)
                {
                    T tmp = base[i];
                    base[i] = base[j];
                    base[j] = tmp;
                    i++;
                    j--;
                }
            } 
            while (i <= j);
            
            if (i < r)
            {
                temp.elem1 = i;
                temp.elem2 = r;
                work_stack.Push(temp);
            }
            r = j;
        }
        while (l < r);
    }
    while (!work_stack.Empty());
}

#endif ///_STACK_H_
