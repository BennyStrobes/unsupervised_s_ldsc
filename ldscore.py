from __future__ import division
import numpy as np
import bitarray as ba
import pdb


def getBlockLefts(coords, max_dist):
    '''
    Converts coordinates + max block length to the a list of coordinates of the leftmost
    SNPs to be included in blocks.

    Parameters
    ----------
    coords : array
        Array of coordinates. Must be sorted.
    max_dist : float
        Maximum distance between SNPs included in the same window.

    Returns
    -------
    block_left : 1D np.ndarray with same length as block_left
        block_left[j] :=  min{k | dist(j, k) < max_dist}.

    '''
    M = len(coords)
    j = 0
    block_left = np.zeros(M)
    for i in xrange(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1

        block_left[i] = j

    return block_left


def block_left_to_right(block_left):
    '''
    Converts block lefts to block rights.

    Parameters
    ----------
    block_left : array
        Array of block lefts.

    Returns
    -------
    block_right : 1D np.ndarray with same length as block_left
        block_right[j] := max {k | block_left[k] <= j}

    '''
    M = len(block_left)
    j = 0
    block_right = np.zeros(M)
    for i in xrange(M):
        while j < M and block_left[j] <= i:
            j += 1

        block_right[i] = j

    return block_right


class __GenotypeArrayInMemory__(object):
    '''
    Parent class for various classes containing interfaces for files with genotype
    matrices, e.g., plink .bed files, etc
    '''
    def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None):
        self.m = len(snp_list.IDList)
        self.n = n
        self.keep_snps = keep_snps
        self.keep_indivs = keep_indivs
        self.df = np.array(snp_list.df[['CHR', 'SNP', 'BP', 'CM']])
        self.colnames = ['CHR', 'SNP', 'BP', 'CM']
        self.mafMin = mafMin if mafMin is not None else 0
        self._currentSNP = 0
        (self.nru, self.geno) = self.__read__(fname, self.m, n)
        # filter individuals
        if keep_indivs is not None:
            keep_indivs = np.array(keep_indivs, dtype='int')
            if np.any(keep_indivs > self.n):
                raise ValueError('keep_indivs indices out of bounds')

            (self.geno, self.m, self.n) = self.__filter_indivs__(self.geno, keep_indivs, self.m,
                self.n)

            if self.n > 0:
                print 'After filtering, {n} individuals remain'.format(n=self.n)
            else:
                raise ValueError('After filtering, no individuals remain')

        # filter SNPs
        if keep_snps is not None:
            keep_snps = np.array(keep_snps, dtype='int')
            if np.any(keep_snps > self.m):  # if keep_snps is None, this returns False
                raise ValueError('keep_snps indices out of bounds')

        (self.geno, self.m, self.n, self.kept_snps, self.freq) = self.__filter_snps_maf__(
            self.geno, self.m, self.n, self.mafMin, keep_snps)

        if self.m > 0:
            print 'After filtering, {m} SNPs remain'.format(m=self.m)
        else:
            raise ValueError('After filtering, no SNPs remain')

        self.df = self.df[self.kept_snps, :]
        self.maf = np.minimum(self.freq, np.ones(self.m)-self.freq)
        self.sqrtpq = np.sqrt(self.freq*(np.ones(self.m)-self.freq))
        self.df = np.c_[self.df, self.maf]
        self.colnames.append('MAF')

    def __read__(self, fname, m, n):
        raise NotImplementedError

    def __filter_indivs__(geno, keep_indivs, m, n):
        raise NotImplementedError

    def __filter_maf_(geno, m, n, maf):
        raise NotImplementedError

    def ldScoreVarBlocks(self, block_left, c, annot=None):
        '''Computes an unbiased estimate of L2(j) for j=1,..,M.'''
        func = lambda x: self.__l2_unbiased__(x, self.n)
        snp_getter = self.nextSNPs
        return self.__corSumVarBlocks__(block_left, c, func, snp_getter, annot)

    def generate_pairwise_lds(self, block_left, snp_names_matrix, pairwise_ld_score_output_file):
        '''Computes an unbiased estimate of L2(j) for j=1,..,M.'''
        func = lambda x: self.__l2_unbiased__(x, self.n)
        c = 1
        self._currentSNP = 0
        snp_getter = self.nextSNPs
        self.__pairwiseCorBlocks__(block_left, c, func, snp_getter, snp_names_matrix, pairwise_ld_score_output_file)

    def generate_pairwise_lds_slow(self, block_left, block_right):
        func = lambda x: self.__l2_unbiased__(x, self.n)
        ldscs = []
        pairwise_lds = []
        pairwise_ld_snps = []

        for snp_index in range(self.m):
            snp_block_left = block_left[snp_index]
            snp_block_right = block_right[snp_index]
            neighbor_snp_genotype = self.snp_slice(snp_block_left, snp_block_right)
            snp_genotype = self.snp_slice(snp_index, snp_index+1)
            corz = np.dot(snp_genotype.T, neighbor_snp_genotype / self.n)
            corz = func(corz)
            ldscs.append(np.sum(corz))
            pairwise_lds.append(corz[0,:])
            pairwise_ld_snps.append(np.arange(snp_block_left, snp_block_right))
        return np.asarray(ldscs), pairwise_lds, pairwise_ld_snps

    def ldScoreBlockJackknife(self, block_left, c, annot=None, jN=10):
        func = lambda x: np.square(x)
        snp_getter = self.nextSNPs
        return self.__corSumBlockJackknife__(block_left, c, func, snp_getter, annot, jN)

    def __l2_unbiased__(self, x, n):
        denom = n-2 if n > 2 else n  # allow n<2 for testing purposes
        sq = np.square(x)
        return sq - (1-sq) / denom


    # general methods for calculating sums of Pearson correlation coefficients
    def __pairwiseCorBlocks__(self, block_left, c, func, snp_getter, snp_names_matrix, pairwise_ld_score_output_file):
        '''
        Parameters
        ----------
        block_left : np.ndarray with shape (M, )
            block_left[i] = index of leftmost SNP included in LD Score of SNP i.
            if c > 1, then only entries that are multiples of c are examined, and it is
            assumed that block_left[a*c+i] = block_left[a*c], except at
            the beginning of the chromosome where the 0th SNP is included in the window.

        c : int
            Chunk size.
        func : function
            Function to be applied to the genotype correlation matrix. Before dotting with
            annot. Examples: for biased L2, np.square. For biased L4,
            lambda x: np.square(np.square(x)). For L1, lambda x: x.
        snp_getter : function(int)
            The method to be used to get the next SNPs (normalized genotypes? Normalized
            genotypes with the minor allele as reference allele? etc)
        annot: numpy array with shape (m,n_a)
            SNP annotations.

        Returns
        -------
        cor_sum : np.ndarray with shape (M, num_annots)
            Estimates.

        '''
        m, n = self.m, self.n
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / c)*c
        #annot = np.ones((m, 1))
         
        num_snps = snp_names_matrix.shape[0]
        t = open(pairwise_ld_score_output_file,'w')
        t.write('CHR\tSNP\tBP\tL2\tpairwise_snp_indices\tpairwise_lds\n')

        #n_a = annot.shape[1]  # number of annotations
        cor_sum = np.zeros(m)
        pairwise_lds = []
        pairwise_ld_snps = []
        for m_iter in range(m):
            pairwise_lds.append([])
            pairwise_ld_snps.append([])
        # b = index of first SNP for which SNP 0 is not included in LD Score
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b/c)*c)  # round up to a multiple of c
        if b > m:
            c = 1
            b = m
        l_A = 0  # l_A := index of leftmost SNP in matrix A
        A = snp_getter(b)
        rfuncAB = np.zeros((b, c))
        rfuncBB = np.zeros((c, c))
        # chunk inside of block
        for l_B in xrange(0, b, c):  # l_B := index of leftmost SNP in matrix B
            B = A[:, l_B:l_B+c]
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)

            cor_sum[l_A:l_A+b] += rfuncAB[:,0]
            pairwise_lds[l_B].append(np.copy(rfuncAB[:,0]))
            pairwise_ld_snps[l_B].append(np.arange(b))
        # chunk to right of block
        b0 = b
        md = int(c*np.floor(m/c))
        end = md + 1 if md != m else md
        for l_B in xrange(b0, end, c):
            # check if the annot matrix is all zeros for this block + chunk
            # this happens w/ sparse categories (i.e., pathways)
            # update the block
            old_b = b
            b = int(block_sizes[l_B])
            if l_B > b0 and b > 0:
                # block_size can't increase more than c
                # block_size can't be less than c unless it is zero
                # both of these things make sense
                A = np.hstack((A[:, old_b-b+c:old_b], B))
                l_A += old_b-b+c
            elif l_B == b0 and b > 0:
                A = A[:, b0-b:b0]
                l_A = b0-b
            elif b == 0:  # no SNPs to left in window, e.g., after a sequence gap
                A = np.array(()).reshape((n, 0))
                l_A = l_B
            if l_B == md:
                c = m - md
                rfuncAB = np.zeros((b, c))
                rfuncBB = np.zeros((c, c))
            if b != old_b:
                rfuncAB = np.zeros((b, c))

            B = snp_getter(c)
            #p1 = np.all(annot[l_A:l_A+b, :] == 0)
            #p2 = np.all(annot[l_B:l_B+c, :] == 0)
            #if p1 and p2:
            #    continue

            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            if b != 0:
                cor_sum[l_A:l_A+b] += rfuncAB[:,0]
                counter = 0
                for index in range(l_A, (l_A+b)):
                    if index > l_B:
                        print('assumption error')
                        pdb.set_trace()
                    if index < block_left[l_B]:
                        print('assumption error')
                        pdb.set_trace()
                    pairwise_lds[index].append(rfuncAB[counter,0])
                    pairwise_ld_snps[index].append(l_B)
                    counter = counter + 1
                cor_sum[l_B] += np.sum(rfuncAB)
                pairwise_lds[l_B].append(np.copy(rfuncAB[:,0]))
                pairwise_ld_snps[l_B].append(range(l_A, l_A+b))
            np.dot(B.T, B / n, out=rfuncBB)
            rfuncBB = func(rfuncBB)
            cor_sum[l_B] += rfuncBB[0,0]
            pairwise_lds[l_B].append(rfuncBB[0,0])
            pairwise_ld_snps[l_B].append(l_B)

            if np.mod(l_B, 5000) == 0:
                for snp_num in range(num_snps):
                    if snp_num >= block_left[l_B] - 2:
                        continue
                    if pairwise_ld_snps[snp_num] == []:
                        continue
                    # Print this snp to output
                    t.write(str(snp_names_matrix[snp_num,0]) + '\t' + snp_names_matrix[snp_num,1] + '\t' + str(snp_names_matrix[snp_num,2]) + '\t')
                    t.write(str(cor_sum[snp_num]) + '\t')
                    t.write(','.join(np.hstack(pairwise_ld_snps[snp_num]).astype(str)) + '\t')
                    t.write(','.join(np.hstack(pairwise_lds[snp_num]).astype(str)) + '\n')
                    pairwise_ld_snps[snp_num] = []
                    pairwise_lds[snp_num] = []
        for snp_num in range(num_snps):
            if pairwise_ld_snps[snp_num] == []:
                continue
            # Print this snp to output
            t.write(str(snp_names_matrix[snp_num,0]) + '\t' + snp_names_matrix[snp_num,1] + '\t' + str(snp_names_matrix[snp_num,2]) + '\t')
            t.write(str(cor_sum[snp_num]) + '\t')
            t.write(','.join(np.hstack(pairwise_ld_snps[snp_num]).astype(str)) + '\t')
            t.write(','.join(np.hstack(pairwise_lds[snp_num]).astype(str)) + '\n')
            pairwise_ld_snps[snp_num] = []
            pairwise_lds[snp_num] = []
        t.close()


        '''
        # Clean up data
        cor_sum2 = []
        for snp_index in range(self.m):
            pairwise_lds[snp_index] = np.hstack(pairwise_lds[snp_index])
            pairwise_ld_snps[snp_index] = np.hstack(pairwise_ld_snps[snp_index])
            cor_sum2.append(np.sum(pairwise_lds[snp_index]))
        cor_sum = np.asarray(cor_sum)
        cor_sum2 = np.asarray(cor_sum2)
        if np.max(cor_sum - cor_sum2) > 1e-5:
            print('assumption error')
            pdb.set_trace()
        if np.abs(np.min(cor_sum-cor_sum2)) > 1e-5:
            print('assumption error')
            pdb.set_trace()
        if snp_names_matrix.shape[0] != len(pairwise_ld_snps):
            print('assumption eror')
            pdb.set_trace()
        num_snps = snp_names_matrix.shape[0]
        t = open(pairwise_ld_score_output_file,'w')
        t.write('CHR\tSNP\tBP\tL2\tpairwise_snp_indices\tpairwise_lds\n')
        for snp_num in range(num_snps):
            t.write(str(snp_names_matrix[snp_num,0]) + '\t' + snp_names_matrix[snp_num,1] + '\t' + str(snp_names_matrix[snp_num,2]) + '\t')
            t.write(str(cor_sum[snp_num]) + '\t')
            t.write(','.join(pairwise_ld_snps[snp_num].astype(str)) + '\t')
            t.write(','.join(pairwise_lds[snp_num].astype(str)) + '\n')
        t.close()
        #return cor_sum, pairwise_lds, pairwise_ld_snps
        '''

    # general methods for calculating sums of Pearson correlation coefficients
    def __pairwiseCorBlocks_memory_hog__(self, block_left, c, func, snp_getter, snp_names_matrix, pairwise_ld_score_output_file):
        '''
        Parameters
        ----------
        block_left : np.ndarray with shape (M, )
            block_left[i] = index of leftmost SNP included in LD Score of SNP i.
            if c > 1, then only entries that are multiples of c are examined, and it is
            assumed that block_left[a*c+i] = block_left[a*c], except at
            the beginning of the chromosome where the 0th SNP is included in the window.

        c : int
            Chunk size.
        func : function
            Function to be applied to the genotype correlation matrix. Before dotting with
            annot. Examples: for biased L2, np.square. For biased L4,
            lambda x: np.square(np.square(x)). For L1, lambda x: x.
        snp_getter : function(int)
            The method to be used to get the next SNPs (normalized genotypes? Normalized
            genotypes with the minor allele as reference allele? etc)
        annot: numpy array with shape (m,n_a)
            SNP annotations.

        Returns
        -------
        cor_sum : np.ndarray with shape (M, num_annots)
            Estimates.

        '''
        m, n = self.m, self.n
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / c)*c
        #annot = np.ones((m, 1))


        #n_a = annot.shape[1]  # number of annotations
        cor_sum = np.zeros(m)
        pairwise_lds = []
        pairwise_ld_snps = []
        for m_iter in range(m):
            pairwise_lds.append([])
            pairwise_ld_snps.append([])
        # b = index of first SNP for which SNP 0 is not included in LD Score
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b/c)*c)  # round up to a multiple of c
        if b > m:
            c = 1
            b = m
        l_A = 0  # l_A := index of leftmost SNP in matrix A
        A = snp_getter(b)
        rfuncAB = np.zeros((b, c))
        rfuncBB = np.zeros((c, c))
        # chunk inside of block
        for l_B in xrange(0, b, c):  # l_B := index of leftmost SNP in matrix B
            B = A[:, l_B:l_B+c]
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)

            cor_sum[l_A:l_A+b] += rfuncAB[:,0]
            pairwise_lds[l_B].append(np.copy(rfuncAB[:,0]))
            pairwise_ld_snps[l_B].append(np.arange(b))
        # chunk to right of block
        b0 = b
        md = int(c*np.floor(m/c))
        end = md + 1 if md != m else md
        for l_B in xrange(b0, end, c):
            print(l_B)
            # check if the annot matrix is all zeros for this block + chunk
            # this happens w/ sparse categories (i.e., pathways)
            # update the block
            old_b = b
            b = int(block_sizes[l_B])
            if l_B > b0 and b > 0:
                # block_size can't increase more than c
                # block_size can't be less than c unless it is zero
                # both of these things make sense
                A = np.hstack((A[:, old_b-b+c:old_b], B))
                l_A += old_b-b+c
            elif l_B == b0 and b > 0:
                A = A[:, b0-b:b0]
                l_A = b0-b
            elif b == 0:  # no SNPs to left in window, e.g., after a sequence gap
                A = np.array(()).reshape((n, 0))
                l_A = l_B
            if l_B == md:
                c = m - md
                rfuncAB = np.zeros((b, c))
                rfuncBB = np.zeros((c, c))
            if b != old_b:
                rfuncAB = np.zeros((b, c))

            B = snp_getter(c)
            #p1 = np.all(annot[l_A:l_A+b, :] == 0)
            #p2 = np.all(annot[l_B:l_B+c, :] == 0)
            #if p1 and p2:
            #    continue

            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            if b != 0:
                cor_sum[l_A:l_A+b] += rfuncAB[:,0]
                counter = 0
                for index in range(l_A, (l_A+b)):
                    if index > l_B:
                        print('assumption error')
                        pdb.set_trace()
                    if index < block_left[l_B]:
                        print('assumption error')
                        pdb.set_trace()
                    pairwise_lds[index].append(rfuncAB[counter,0])
                    pairwise_ld_snps[index].append(l_B)
                    counter = counter + 1
                cor_sum[l_B] += np.sum(rfuncAB)
                pairwise_lds[l_B].append(np.copy(rfuncAB[:,0]))
                pairwise_ld_snps[l_B].append(range(l_A, l_A+b))
            np.dot(B.T, B / n, out=rfuncBB)
            rfuncBB = func(rfuncBB)
            cor_sum[l_B] += rfuncBB[0,0]
            pairwise_lds[l_B].append(rfuncBB[0,0])
            pairwise_ld_snps[l_B].append(l_B)

        # Clean up data
        cor_sum2 = []
        for snp_index in range(self.m):
            pairwise_lds[snp_index] = np.hstack(pairwise_lds[snp_index])
            pairwise_ld_snps[snp_index] = np.hstack(pairwise_ld_snps[snp_index])
            cor_sum2.append(np.sum(pairwise_lds[snp_index]))
        cor_sum = np.asarray(cor_sum)
        cor_sum2 = np.asarray(cor_sum2)
        if np.max(cor_sum - cor_sum2) > 1e-5:
            print('assumption error')
            pdb.set_trace()
        if np.abs(np.min(cor_sum-cor_sum2)) > 1e-5:
            print('assumption error')
            pdb.set_trace()
        if snp_names_matrix.shape[0] != len(pairwise_ld_snps):
            print('assumption eror')
            pdb.set_trace()
        num_snps = snp_names_matrix.shape[0]
        t = open(pairwise_ld_score_output_file,'w')
        t.write('CHR\tSNP\tBP\tL2\tpairwise_snp_indices\tpairwise_lds\n')
        for snp_num in range(num_snps):
            t.write(str(snp_names_matrix[snp_num,0]) + '\t' + snp_names_matrix[snp_num,1] + '\t' + str(snp_names_matrix[snp_num,2]) + '\t')
            t.write(str(cor_sum[snp_num]) + '\t')
            t.write(','.join(pairwise_ld_snps[snp_num].astype(str)) + '\t')
            t.write(','.join(pairwise_lds[snp_num].astype(str)) + '\n')
        t.close()
        #return cor_sum, pairwise_lds, pairwise_ld_snps



    # general methods for calculating sums of Pearson correlation coefficients
    def __corSumVarBlocks__(self, block_left, c, func, snp_getter, annot=None):
        '''
        Parameters
        ----------
        block_left : np.ndarray with shape (M, )
            block_left[i] = index of leftmost SNP included in LD Score of SNP i.
            if c > 1, then only entries that are multiples of c are examined, and it is
            assumed that block_left[a*c+i] = block_left[a*c], except at
            the beginning of the chromosome where the 0th SNP is included in the window.

        c : int
            Chunk size.
        func : function
            Function to be applied to the genotype correlation matrix. Before dotting with
            annot. Examples: for biased L2, np.square. For biased L4,
            lambda x: np.square(np.square(x)). For L1, lambda x: x.
        snp_getter : function(int)
            The method to be used to get the next SNPs (normalized genotypes? Normalized
            genotypes with the minor allele as reference allele? etc)
        annot: numpy array with shape (m,n_a)
            SNP annotations.

        Returns
        -------
        cor_sum : np.ndarray with shape (M, num_annots)
            Estimates.

        '''
        m, n = self.m, self.n
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / c)*c
        if annot is None:
            annot = np.ones((m, 1))
        else:
            annot_m = annot.shape[0]
            if annot_m != self.m:
                raise ValueError('Incorrect number of SNPs in annot')

        n_a = annot.shape[1]  # number of annotations
        cor_sum = np.zeros((m, n_a))
        # b = index of first SNP for which SNP 0 is not included in LD Score
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b/c)*c)  # round up to a multiple of c
        if b > m:
            c = 1
            b = m
        l_A = 0  # l_A := index of leftmost SNP in matrix A
        A = snp_getter(b)
        rfuncAB = np.zeros((b, c))
        rfuncBB = np.zeros((c, c))
        # chunk inside of block
        for l_B in xrange(0, b, c):  # l_B := index of leftmost SNP in matrix B
            B = A[:, l_B:l_B+c]
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A:l_A+b, :] += np.dot(rfuncAB, annot[l_B:l_B+c, :])
        # chunk to right of block
        b0 = b
        md = int(c*np.floor(m/c))
        end = md + 1 if md != m else md
        for l_B in xrange(b0, end, c):
            # check if the annot matrix is all zeros for this block + chunk
            # this happens w/ sparse categories (i.e., pathways)
            # update the block
            old_b = b
            b = int(block_sizes[l_B])
            if l_B > b0 and b > 0:
                # block_size can't increase more than c
                # block_size can't be less than c unless it is zero
                # both of these things make sense
                A = np.hstack((A[:, old_b-b+c:old_b], B))
                l_A += old_b-b+c
            elif l_B == b0 and b > 0:
                A = A[:, b0-b:b0]
                l_A = b0-b
            elif b == 0:  # no SNPs to left in window, e.g., after a sequence gap
                A = np.array(()).reshape((n, 0))
                l_A = l_B
            if l_B == md:
                c = m - md
                rfuncAB = np.zeros((b, c))
                rfuncBB = np.zeros((c, c))
            if b != old_b:
                rfuncAB = np.zeros((b, c))

            B = snp_getter(c)
            p1 = np.all(annot[l_A:l_A+b, :] == 0)
            p2 = np.all(annot[l_B:l_B+c, :] == 0)
            if p1 and p2:
                continue

            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A:l_A+b, :] += np.dot(rfuncAB, annot[l_B:l_B+c, :])
            cor_sum[l_B:l_B+c, :] += np.dot(annot[l_A:l_A+b, :].T, rfuncAB).T
            np.dot(B.T, B / n, out=rfuncBB)
            rfuncBB = func(rfuncBB)
            cor_sum[l_B:l_B+c, :] += np.dot(rfuncBB, annot[l_B:l_B+c, :])

        return cor_sum


class PlinkBEDFile(__GenotypeArrayInMemory__):
    '''
    Interface for Plink .bed format
    '''
    def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None):
        self._bedcode = {
            2: ba.bitarray('11'),
            9: ba.bitarray('10'),
            1: ba.bitarray('01'),
            0: ba.bitarray('00')
            }

        __GenotypeArrayInMemory__.__init__(self, fname, n, snp_list, keep_snps=keep_snps,
            keep_indivs=keep_indivs, mafMin=mafMin)

    def __read__(self, fname, m, n):
        if not fname.endswith('.bed'):
            raise ValueError('.bed filename must end in .bed')

        fh = open(fname, 'rb')
        magicNumber = ba.bitarray(endian="little")
        magicNumber.fromfile(fh, 2)
        bedMode = ba.bitarray(endian="little")
        bedMode.fromfile(fh, 1)
        e = (4 - n % 4) if n % 4 != 0 else 0
        nru = n + e
        self.nru = nru
        # check magic number
        if magicNumber != ba.bitarray('0011011011011000'):
            raise IOError("Magic number from Plink .bed file not recognized")

        if bedMode != ba.bitarray('10000000'):
            raise IOError("Plink .bed file must be in default SNP-major mode")

        # check file length
        self.geno = ba.bitarray(endian="little")
        self.geno.fromfile(fh)
        self.__test_length__(self.geno, self.m, self.nru)
        return (self.nru, self.geno)

    def __test_length__(self, geno, m, nru):
        exp_len = 2*m*nru
        real_len = len(geno)
        if real_len != exp_len:
            s = "Plink .bed file has {n1} bits, expected {n2}"
            raise IOError(s.format(n1=real_len, n2=exp_len))

    def __filter_indivs__(self, geno, keep_indivs, m, n):
        n_new = len(keep_indivs)
        e = (4 - n_new % 4) if n_new % 4 != 0 else 0
        nru_new = n_new + e
        nru = self.nru
        z = ba.bitarray(m*2*nru_new, endian="little")
	z.setall(0)
        for e, i in enumerate(keep_indivs):
            z[2*e::2*nru_new] = geno[2*i::2*nru]
            z[2*e+1::2*nru_new] = geno[2*i+1::2*nru]

        self.nru = nru_new
        return (z, m, n_new)

    def __filter_snps_maf__(self, geno, m, n, mafMin, keep_snps):
        '''
        Credit to Chris Chang and the Plink2 developers for this algorithm
        Modified from plink_filter.c
        https://github.com/chrchang/plink-ng/blob/master/plink_filter.c

        Genotypes are read forwards (since we are cheating and using endian="little")

        A := (genotype) & 1010...
        B := (genotype) & 0101...
        C := (A >> 1) & B

        Then

        a := A.count() = missing ct + hom major ct
        b := B.count() = het ct + hom major ct
        c := C.count() = hom major ct

        Which implies that

        missing ct = a - c
        # of indivs with nonmissing genotype = n - a + c
        major allele ct = b + c
        major allele frequency = (b+c)/(2*(n-a+c))
        het ct + missing ct = a + b - 2*c

        Why does bitarray not have >> ????

        '''
        nru = self.nru
        m_poly = 0
        y = ba.bitarray()
        if keep_snps is None:
            keep_snps = xrange(m)
        kept_snps = []
        freq = []
        for e, j in enumerate(keep_snps):
            z = geno[2*nru*j:2*nru*(j+1)]
            A = z[0::2]
            a = A.count()
            B = z[1::2]
            b = B.count()
            c = (A & B).count()
            major_ct = b + c  # number of copies of the major allele
            n_nomiss = n - a + c  # number of individuals with nonmissing genotypes
            f = major_ct / (2*n_nomiss) if n_nomiss > 0 else 0
            het_miss_ct = a+b-2*c  # remove SNPs that are only either het or missing
            if np.minimum(f, 1-f) > mafMin and het_miss_ct < n:
                freq.append(f)
                y += z
                m_poly += 1
                kept_snps.append(j)

        return (y, m_poly, n, kept_snps, freq)

    def nextSNPs(self, b, minorRef=None):
        '''
        Unpacks the binary array of genotypes and returns an n x b matrix of floats of
        normalized genotypes for the next b SNPs, where n := number of samples.

        Parameters
        ----------
        b : int
            Number of SNPs to return.
        minorRef: bool, default None
            Should we flip reference alleles so that the minor allele is the reference?
            (This is useful for computing l1 w.r.t. minor allele).

        Returns
        -------
        X : np.array with dtype float64 with shape (n, b), where n := number of samples
            Matrix of genotypes normalized to mean zero and variance one. If minorRef is
            not None, then the minor allele will be the positive allele (i.e., two copies
            of the minor allele --> a positive number).

        '''

        try:
            b = int(b)
            if b <= 0:
                raise ValueError("b must be > 0")
        except TypeError:
            raise TypeError("b must be an integer")

        if self._currentSNP + b > self.m:
            s = '{b} SNPs requested, {k} SNPs remain'
            raise ValueError(s.format(b=b, k=(self.m-self._currentSNP)))

        c = self._currentSNP
        n = self.n
        nru = self.nru
        slice = self.geno[2*c*nru:2*(c+b)*nru]
        X = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, nru)).T
        X = X[0:n, :]
        Y = np.zeros(X.shape)
        for j in xrange(0, b):
            newsnp = X[:, j]
            ii = newsnp != 9
            avg = np.mean(newsnp[ii])
            newsnp[np.logical_not(ii)] = avg
            denom = np.std(newsnp)
            if denom == 0:
                denom = 1

            if minorRef is not None and self.freq[self._currentSNP + j] > 0.5:
                denom = denom*-1

            Y[:, j] = (newsnp - avg) / denom

        self._currentSNP += b
        return Y
    def snp_slice(self, start, end, minorRef=None):
        '''
        Unpacks the binary array of genotypes and returns an n x b matrix of floats of
        normalized genotypes for the next b SNPs, where n := number of samples.

        Parameters
        ----------
        b : int
            Number of SNPs to return.
        minorRef: bool, default None
            Should we flip reference alleles so that the minor allele is the reference?
            (This is useful for computing l1 w.r.t. minor allele).

        Returns
        -------
        X : np.array with dtype float64 with shape (n, b), where n := number of samples
            Matrix of genotypes normalized to mean zero and variance one. If minorRef is
            not None, then the minor allele will be the positive allele (i.e., two copies
            of the minor allele --> a positive number).

        '''

        try:
            start = int(start)
            if start < 0:
                raise ValueError("start must be > 0")
        except TypeError:
            raise TypeError("start must be an integer")
        try:
            end = int(end)
            if end <= 0:
                raise ValueError("end must be > 0")
        except TypeError:
            raise TypeError("end must be an integer")

        if start > self.m:
            raise ValueError("start is greater than number of snps")

        if end > self.m:
            raise ValueError("start is greater than number of snps")

        if end < start:
            raise ValueError("start is less than end")

        n = self.n
        nru = self.nru
        b = end - start
        slice = self.geno[2*start*nru:2*(end)*nru]
        X = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, nru)).T
        X = X[0:n, :]
        Y = np.zeros(X.shape)
        for j in xrange(0, b):
            newsnp = X[:, j]
            ii = newsnp != 9
            avg = np.mean(newsnp[ii])
            newsnp[np.logical_not(ii)] = avg
            denom = np.std(newsnp)
            if denom == 0:
                denom = 1

            if minorRef is not None and self.freq[start + j] > 0.5:
                denom = denom*-1

            Y[:, j] = (newsnp - avg) / denom

        return Y
