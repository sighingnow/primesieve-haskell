-----------------------------------------------------------
-- |
-- module:                      Math.Linear.Matrix
-- copyright:                   (c) 2016-2017 HE, Tao
-- license:                     MIT
-- maintainer:                  sighingnow@gmail.com
--
-- Matrices representation in Haskell.
--
{-# OPTIONS_GHC -Wno-orphans #-}
{-# OPTIONS_GHC -Wno-inline-rule-shadowing #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}

module Math.Linear.Matrix where
  -- ( Mat(..)
  --  -- * Constructors.
  -- , zeros
  -- , ones
  -- , identity
  -- , random
  -- , permute
  -- , replicate
  -- , matrix
  --  -- * Convertion between list and matrix.
  -- , diag
  -- , toList
  -- , fromList
  --  -- * Predicates.
  -- , empty
  -- , null
  -- , valid
  -- , square
  --  -- * Data accessor.
  -- , (!)
  -- , diagonal
  -- , rows
  -- , rowAt
  -- , columns
  -- , columnAt
  --  -- * Searching
  -- , elem
  -- , find
  --  -- * Folding
  -- , all
  -- , any
  -- , foldl'
  --  -- * Matrix algbera.
  -- , shift
  -- , times
  -- , add
  -- , minus
  -- , mult
  -- , division
  -- , (.*)
  -- , dot
  -- , dot'
  -- , pow
  -- , safeDot
  -- , inner
  -- , inner'
  --  -- * Matrix reshape operations.
  -- , (==>)
  -- , reshape
  -- , transpose
  -- , lower
  -- , upper
  -- , concatenate
  -- , (<|>)
  -- , (<->)
  --  -- * Matrix properties.
  -- , shape
  -- , sum
  -- , product
  -- , minimum
  -- , maximum
  -- , mean
  --  -- Generic utilities.
  -- , map
  -- , filter
  -- , zipWith
  -- , force
  --  -- Mutable matrix and raw pointers.
  -- , unsafeWith
  -- , unsafeThaw
  -- , unsafeFreeze
  -- , unsafeOp
  -- , unsafeUnaryOp
  -- , unsafeBinaryOp
  -- , unsafeDecomposeOp)
  -- where

import Foundation
import Foundation.Array.Internal (withPtr, withMutablePtr)
import Foundation.Class.Storable
import Foundation.Collection
import Foundation.Primitive

import GHC.Num (Num)
import Foreign.Marshal.Alloc (alloca)
import System.IO.Unsafe (unsafePerformIO)

import qualified Math.Linear.Internal as I
import qualified Math.Linear.Matrix.Mutable as Mutable

data Mat a = M
    { row :: {-# UNPACK #-}!Int -- ^ rows
    , column :: {-# UNPACK #-}!Int -- ^ columns
    , vect :: {-# UNPACK #-}!(UArray a) -- ^ data in plain vector.
    } deriving (Eq, Show)

-- | Construct an empty matrix.
empty :: PrimType a => Mat a
empty = M 0 0 mempty

{-# INLINE empty #-}

-- | Verify matrix dimensions and memory layout
valid :: PrimType a => Mat a -> Bool
valid M{..} = row >= 0 && column >= 0 && length vect == CountOf (row * column)

{-# INLINE valid #-}

-- | If the matrix is a square matrix.
square :: Mat a -> Bool
square M{..} = row == column

{-# INLINE square #-}

-- | Construct a matrix with all zeros.
zeros :: (PrimType a, Num a) => Int -> Int -> Mat a
zeros r c = replicate' r c 0

{-# INLINE zeros #-}

-- | Construct a matrix with all ones.
ones :: (PrimType a, Num a) => Int -> Int -> Mat a
ones r c = replicate' r c 1

{-# INLINE ones #-}

-- | Construct a identity matrix, square is not required.
identity :: I.Elem a
    => Int -> Int -> Mat a
identity r c = unsafePerformIO $ do
    m <- Mutable.zeros r c
    Mutable.unsafeWith m $ \xs r' c' ->
        I.call $ I.identity xs r' c'
    unsafeFreeze m

{-# INLINE identity #-}

-- | Construct a random matrix.
random :: I.Elem a
    => Int -> Int -> IO (Mat a)
random r c = do
    m <- Mutable.new r c
    Mutable.unsafeWith m $ \xs r' c' ->
        I.call' $ I.random_ xs r' c'
    unsafeFreeze m

{-# INLINE random #-}

-- | Permutation matrix, assume that /P/ is the permutation matrix, then,
--
-- * exchange two columns: /x' = x A/.
--
-- * exchange two rows: /x' = A x/.
permute :: I.Elem a
    => Int -- ^ size of matrix, /n x n/.
    -> Int
    -> Int
    -> Mat a
permute n i j = unsafePerformIO $ do
    m <- Mutable.identity n n
    Mutable.unsafeWrite m i i 0
    Mutable.unsafeWrite m j j 0
    Mutable.unsafeWrite m i j 1
    Mutable.unsafeWrite m j i 1
    unsafeFreeze m

{-# INLINE permute #-}

-- | Construct a matrix with given constant.
replicate' :: PrimType a
    => Int -> Int -> a -> Mat a
replicate' r c = M r c . replicate (CountOf (r * c))

{-# INLINE replicate' #-}

-- | Construct matrix with given generate function.
matrix :: PrimType a
    => Int -> Int -> (Int -> Int -> a) -> Mat a
matrix r c func = fromList' r c [ func i j
                                | j <- [0 .. c - 1]
                                , i <- [0 .. r - 1] ]

{-# INLINE matrix #-}

-- | Construct a diagonal matrix from given vector, elements in vector will be diagonal elements of the matrix.
diag :: I.Elem a
    => UArray a -> Mat a
diag xs = unsafePerformIO $ do
    let nlen = fromCount (length xs)
    m <- Mutable.new nlen nlen
    Mutable.unsafeWith m $ \vect r c ->
        withPtr xs $ \p ->
            I.call $ I.diag vect r c p
    unsafeFreeze m

{-# INLINE [1] diag #-}

-- | Construct matrix from given list and size, assume that /row * column == length xs/.
fromList' :: PrimType a
    => Int -- ^ rows
    -> Int -- ^ columns
    -> [a] -- ^ values
    -> Mat a
fromList' r c = M r c . fromListN (r * c)

{-# INLINE [1] fromList' #-}

{-# RULES
"fromlist'/tolist" forall a b m. fromList' a b (toList m) = m
 #-}

{-# RULES
"tolist/fromlist'" forall a b l. toList (fromList' a b l) = l
 #-} -- assume /a * b == length l/

-- | If the vector contains a specific element that statisfy some predicate.
find' :: PrimType a => (a -> Bool) -> Mat a -> Maybe a
find' predicate = find predicate . vect

{-# INLINE find' #-}

-- -- | Left fold strictly, /O(n)/.
-- foldl'
--     :: Storable b
--     => (a -> b -> a) -> a -> Mat b -> a
-- foldl' f x = V.foldl' f x . vect

-- {-# INLINE foldl' #-}

-- | Get specified element from matrix.
at :: PrimType a
    => Mat a -> (Int, Int) -> Maybe a
at M{..} (i, j) = vect ! (Offset (i * column + j))

{-# INLINE at #-}

-- | Unsafely get specified element from matrix.
at' :: PrimType a
    => Mat a -> (Int, Int) -> a
at' M{..} (i, j) = case vect ! (Offset (i * column + j)) of
                     Just v -> v
                     Nothing -> error "Matrix.at': no such element."

{-# INLINE at' #-}

-- | Get elements in diagonal positions, square matrix is not necessary.
diagonal :: I.Elem a
    => Mat a -> UArray a
diagonal m@M{..} = unsafePerformIO $ do
    v <- mutNew (CountOf (min row column))
    unsafeWith m $ \xs r c ->
        withMutablePtr v $ \p ->
            I.call $ I.diagonal p xs r c
    unsafeFreeze v

{-# INLINE [1] diagonal #-}

{-# RULES
"diagonal/diag" forall v. diagonal (diag v) = v
 #-}

-- | Get deminsion of the matrix.
shape :: Mat a -> (Int, Int)
shape M{..} = (row, column)

{-# INLINE [1] shape #-}

-- | Reshape a matrix.
reshape :: Mat a -> (Int, Int) -> Mat a
reshape M{..} (r, c)
    | row * column == r * c = M r c vect
    | otherwise = error "Matrix.reshape: total size of new matrix must be unchanged."

{-# INLINE [1] reshape #-}

{-# RULES
"reshape/reshape" forall a b m. reshape (reshape m a) b =
                  reshape m b
"reshape/shape" forall m. reshape m (shape m) = m
 #-}

-- | Infix operator of reshape.
(==>) :: Mat a -> (Int, Int) -> Mat a
(==>) = reshape

infix 8 ==>

{-# INLINE (==>) #-}

-- -- | Lines of a matrix.
-- rows :: Storable a => Mat a -> [V.Vector a]
-- rows m@M {..} = Prelude.map (rowAt m) [0 .. row - 1]

-- {-# INLINE rows #-}

-- -- | Get specific row.
-- rowAt
--     :: Storable a
--     => Mat a -> Int -> V.Vector a
-- rowAt M {..} r = V.slice (r * column) column vect

-- {-# INLINE rowAt #-}

-- -- | Columns of a matrix.
-- columns
--     :: Storable a
--     => Mat a -> [V.Vector a]
-- columns m@M {..} = Prelude.map (columnAt m) [0 .. column - 1]

-- {-# INLINE columns #-}

-- -- | Get specific column.
-- columnAt
--     :: Storable a
--     => Mat a -> Int -> V.Vector a
-- columnAt M {..} c = V.generate row $ \i -> vect V.! (i * column + c)

-- {-# INLINE columnAt #-}

-- -- | Map an unary function over a matrix.
-- map
--     :: (Storable a, Storable b)
--     => (a -> b) -> Mat a -> Mat b
-- map f M {..} = M row column $ V.map f vect

-- {-# INLINE [1] map #-}

-- {-# RULES
-- "map/map" forall f g m. map f (map g m) = map (f . g) m
--  #-}

-- -- | Filter elements from a matrix ans accumulate these elements into a vector.
-- filter
--     :: Storable a
--     => (a -> Bool) -> Mat a -> V.Vector a
-- filter f = V.filter f . vect

-- {-# INLINE [1] filter #-}

-- {-# RULES
-- "filter/map" forall f g m. filter f (map g m) = filter (f . g) m
--  #-}

-- | Matrix transpose.
transpose :: I.Elem a => Mat a -> Mat a
transpose m@M{..} = unsafeUnaryOp column row I.transpose m -- matrix column row $ \i j -> m ! (j, i)

{-# INLINE [1] transpose #-}

{-# RULES
"transpose/transpose" forall m. transpose (transpose m) = m
 #-}

-- | Lower triangularize the given matrix strictly.
lower :: I.Elem a => Mat a -> Mat a
lower m@M{..} = unsafeUnaryOp row column I.lower m

{-# INLINE lower #-}

-- | Upper triangularize the given matrix strictly.
upper :: I.Elem a => Mat a -> Mat a
upper m@M{..} = unsafeUnaryOp row column I.upper m

{-# INLINE upper #-}

-- -- | Concatenate mutiple matrices into one with row as axis, columns of these operands must be equal.
-- --
-- -- > concatenate [[1,2],[3,4]] [[5,6],[7,8]]
-- -- >    1 2
-- -- >    3 4
-- -- >    5 6
-- -- >    7 8
-- concatenate
--     :: Storable a
--     => [Mat a] -> Mat a
-- concatenate [] = error "Matrix.concatenate: no operand found."
-- concatenate ms
--     | length csc <= 1 = M rs (List.head csc) vects
--     | otherwise = error "Matrix.concatenate: all matrice must have same column size."
--   where
--     rs = Prelude.sum . Prelude.map row $ ms
--     csc = List.concat . List.group . Prelude.map column $ ms
--     vects = V.concat . Prelude.map vect $ ms

-- {-# INLINE concatenate #-}

-- -- | Horizontally join two matrices.
-- --
-- -- > ( A ) <|> ( B ) = ( A | B )
-- --
-- -- Where both matrices /A/ and /B/ have the same number of rows, but no bounds check will be performed.
-- (<|>)
--     :: Storable a
--     => Mat a -> Mat a -> Mat a
-- m1 <|> m2 =
--     matrix (row m1) (c + column m2) $
--     \i j ->
--          if j <= c
--              then m1 ! (i, j)
--              else m2 ! (i, j - c)
--   where
--     c = column m1

-- -- {-# INLINE (<|>) #-}

-- -- | Vertically join two matrices.
-- --
-- -- >                   ( A )
-- -- > ( A ) <-> ( B ) = ( - )
-- -- >                   ( B )
-- --
-- -- Where both matrices /A/ and /B/ have the same number of columns, but no bounds check will be performed.
-- (<->)
--     :: Storable a
--     => Mat a -> Mat a -> Mat a
-- m1 <-> m2 =
--     matrix (r + row m2) (column m1) $
--     \i j ->
--          if i <= r
--              then m1 ! (i, j)
--              else m2 ! (i - r, j)
--   where
--     r = row m1

-- {-# INLINE (<->) #-}

-- -- | Zip two matrix with given function and truncate extra elements both row and column.
-- zipWith
--     :: (Storable a, Storable b, Storable c)
--     => (a -> b -> c) -> Mat a -> Mat b -> Mat c
-- zipWith f m1@(M r1 c1 _) m2@(M r2 c2 _)
--     | r1 == r1 && c1 == c2 = M r1 c1 $ V.zipWith f (vect m1) (vect m2)
--     | r1 <= r2 && c1 <= c2 = matrix r1 c1 $ \i j -> f (m1 ! (i, j)) (m2 ! (i, j))
--     | r1 >= r2 && c1 >= c2 = matrix r2 c2 $ \i j -> f (m1 ! (i, j)) (m2 ! (i, j))
--     | otherwise = error "Mat.zipWith: uncompatible shape of two operands."

-- | The sum of all coefficients of the matrix
sum :: I.Elem a => Mat a -> a
sum = unsafeOp I.sum -- V.sum . vect

{-# INLINE sum #-}

-- | The product of all coefficients of the matrix
product :: I.Elem a => Mat a -> a
product = unsafeOp I.product -- V.product . vect

{-# INLINE product #-}

-- | The mean of all coefficients of the matrix
mean :: I.Elem a => Mat a -> a
mean = unsafeOp I.mean

{-# INLINE mean #-}

-- | Add scalar value to every elements in matrix.
shift :: (I.Elem a, Additive a) => a -> Mat a -> Mat a
shift x m@M {..} = unsafePerformIO $ do
    m' <- Mutable.new row column
    Mutable.unsafeWith m' $ \v' _ _ ->
        unsafeWith m $ \v _ _ ->
            alloca $ \p -> do
                poke p x
                I.call $ I.shift v' p v row column
    unsafeFreeze m'

{-# INLINE [1] shift #-}

{-# RULES
"shift/shift" forall a b m. shift a (shift b m) = shift (a + b) m
 #-}

-- | Multiply scalar value to every elements in matrix.
times :: (I.Elem a, Multiplicative a) => a -> Mat a -> Mat a
times x m@M {..} = unsafePerformIO $ do
    m' <- Mutable.new row column
    Mutable.unsafeWith m' $ \v' _ _ ->
        unsafeWith m $ \v _ _ ->
            alloca $ \p -> do
                poke p x
                I.call $ I.times v' p v row column
    unsafeFreeze m'

{-# INLINE [1] times #-}

{-# RULES
"times/times" forall a b m. times a (times b m) = times (a * b) m
 #-}

-- | Elementwise addition.
add :: I.Elem a => Mat a -> Mat a -> Mat a
add m1 m2 = unsafeBinaryOp (row m1) (column m2) I.add m1 m2

{-# INLINE add #-}

-- | Elementwise substraction.
minus :: I.Elem a => Mat a -> Mat a -> Mat a
minus m1 m2 = unsafeBinaryOp (row m1) (column m2) I.minus m1 m2

{-# INLINE minus #-}

-- | Elementwise multiplication.
mult :: I.Elem a => Mat a -> Mat a -> Mat a
mult m1 m2 = unsafeBinaryOp (row m1) (column m2) I.mult m1 m2

{-# INLINE mult #-}

-- | Elementwise division for real fractions.
division :: I.Elem a => Mat a -> Mat a -> Mat a
division m1 m2 = unsafeBinaryOp (row m1) (column m2) I.division m1 m2

{-# INLINE division #-}

-- | Matrix multiplication WITHOUT bounds check.
dot :: I.Elem a => Mat a -> Mat a -> Mat a
dot m1 m2 = unsafeBinaryOp (row m1) (column m2) I.dot m1 m2

{-# INLINE dot #-}

-- | Infix operator of dot.
(.*) :: I.Elem a => Mat a -> Mat a -> Mat a
(.*) = dot

infix 8 .*

{-# INLINE (.*) #-}

-- | Matrix multiplication with bounds check.
dot' :: I.Elem a
    => Mat a -> Mat a -> Mat a
dot' m1 m2
    | column m1 == row m2 = dot m1 m2
    | otherwise = error "Mat.dot: unmatched shape of two operands."

{-# INLINE dot' #-}

-- -- | Matrix multiplication in pure Haskell without FFI.
-- dot'' :: I.Elem a => Mat a -> Mat a -> Mat a
-- dot'' m1 m2 = matrix (row m1) (column m2) $ \i j -> V.sum $ V.zipWith (*) (rowAt m1 i) (columnAt m2 j)

-- {-# INLINE dot'' #-}

-- | Matrix power with matrix multiplication (dot) as multiply operator, the given matrix must
-- be a square matrix.
pow :: (I.Elem a, Integral b, Ord b, Num b, IDivisible b)
    => Mat a -> b -> Mat a
pow m@M {..} k
    | row /= column = error "Matrix.pow: a square matrix is needed."
    | k < 0 = error "Matrix.pow: negative exponent is not allowed."
    | otherwise = go (diag $ replicate (CountOf row) 1) m k
  where
    go r _ 0 = r
    go r a n = go r' a' n'
      where
        r' = if n `mod` 2 == 0
                then r
                else r `dot` a
        a' = a `dot` a
        n' = n `div` 2

{-# INLINE pow #-}

-- | inner product of two specific columns (column i and column j) WITHOUT bounds check.
inner :: I.Elem a
    => Mat a -> Int -> Int -> a
inner m i j = unsafePerformIO $
    unsafeWith m $ \xs r c ->
        alloca $ \p -> do
            I.call $ I.inner p r xs c i xs c j
            peek p

{-# INLINE inner #-}

-- | inner product of two vectors with same length, but no bounds check will be performed.
inner' :: I.Elem a
    => UArray a -> UArray a -> a
inner' v1 v2 = unsafePerformIO $
    withPtr v1 $ \xs1 ->
        withPtr v2 $ \xs2 ->
            alloca $ \p -> do
                I.call $ I.inner p n xs1 1 0 xs2 1 0
                peek p
  where
    n = fromCount $ min (length v1) (length v2)

{-# INLINE inner' #-}

det :: I.Elem a => Mat a -> a
det m = unsafePerformIO $ do
    unsafeWith m $ \xs r c ->
      alloca $ \res -> do
        I.call $ I.det res xs r c
        peek res

{-# LANGUAGE det #-}

-- | Copy the matrix data and drop extra memory.
force :: PrimType a => Mat a -> Mat a
force m@M{..} = matrix row column $ curry (m `at'`)

{-# INLINE force #-}

-- | Pass a pointer to the matrix's data to the IO action. The data may not be modified through the pointer.
unsafeWith :: (PrimMonad m, I.Elem a)
    => Mat a -> (Ptr a -> Int -> Int -> m c) -> m c
unsafeWith M{..} f = withPtr vect $ \p -> f p row column

{-# INLINE unsafeWith #-}

-- | Unsafe convert a mutable matrix to an immutable one without copying.
-- The mutable matrix may not be used after this operation.
unsafeFreeze' :: (PrimMonad m, I.Elem a)
    => Mutable.MMat a (PrimState m) -> m (Mat a)
unsafeFreeze' Mutable.M{..} = M row column <$> unsafeFreeze vect

{-# INLINE unsafeFreeze' #-}

-- | Unsafely convert an immutable matrix to a mutable one without copying.
-- The immutable matrix may not be used after this operation.
unsafeThaw' :: (PrimMonad m, I.Elem a)
    => Mat a -> m (Mutable.MMat a (PrimState m))
unsafeThaw' M{..} = Mutable.M row column <$> unsafeThaw vect

{-# INLINE unsafeThaw' #-}

-- | Take one matrix as argument and only return a scalar, like sum, mean, max and min.
unsafeOp :: (I.Elem a, I.Elem b)
    => (Ptr b -- result ptr, a scalar value.
        -> Ptr a -- matrix operand
        -> Int -- rows of matrix operand
        -> Int -- columns of matrix operand
        -> Int -- CFFI call status code, non-zero means error
      ) -- ^ op function
    -> Mat a -- ^ operand
    -> b
unsafeOp f m = unsafePerformIO $
    unsafeWith m $ \xs r c ->
        alloca $ \p -> do
            I.call $ f p xs r c
            peek p

{-# INLINE unsafeOp #-}

-- | Take one matrix as it's argument and return a matrix as result, like transpose and inverse.
unsafeUnaryOp :: I.Elem a
    => Int -- ^ target rows
    -> Int -- ^ target columns
    -> (Ptr a -- ^ result ptr, a matrix
        -> Ptr a -- ^ matrix operand ptr
        -> Int -- ^ rows of matrix operand
        -> Int -- ^ columns of matrix operand
        -> Int -- ^ CFFI call status code, non-zero means error
      ) -- ^ op function
    -> Mat a -- ^ matrix operand
    -> Mat a
unsafeUnaryOp r c f m = unsafePerformIO $ do
    m' <- Mutable.zeros r c
    Mutable.unsafeWith m' $ \xs' _ _ -> unsafeWith m $ \xs u v -> I.call $ f xs' xs u v
    unsafeFreeze m'

{-# INLINE unsafeUnaryOp #-}

-- | Take two matrices as arguments and return a matrix as result, like dot, add, multiplication.
unsafeBinaryOp :: I.Elem a
    => Int -- ^ target rows
    -> Int -- ^ target columns
    -> (Ptr a -- ^ result ptr, a matrix, with shape /m x n/
        -> Int -- ^ m
        -> Int -- ^ n
        -> Int -- ^ k
        -> Ptr a -- ^ left matrix operand, with shape /m x k/
        -> Ptr a -- ^ right matrix operand, with shape /k x n/
        -> Int -- ^ CFFI call status code, non-zero means error
      ) -- ^ op function
    -> Mat a -- ^ left operand
    -> Mat a -- ^ right operand
    -> Mat a
unsafeBinaryOp r c f m1 m2 = unsafePerformIO $ do
  m0 <- Mutable.new r c
  Mutable.unsafeWith m0 $
    \vect0 m n ->
      unsafeWith m1 $
        \vect1 _ k -> unsafeWith m2 $ \vect2 _ _ -> I.call $ f vect0 m n k vect1 vect2
  unsafeFreeze' m0

{-# INLINE unsafeBinaryOp #-}

-- | Take one matrix as argument and three matrices as result, like LU decomposition, QR decomposition and SVD.
unsafeDecomposeOp
    :: I.Elem a
    => (Int, Int) -- ^ size of result 1
    -> (Int, Int) -- ^ size of result 2
    -> (Int, Int) -- ^ size of result 3
    -> (Ptr a -- ^ result ptr 1, a matrix
        -> Int -> Int -- ^ size of result 1
        -> Ptr a -- ^ result ptr 2, a matrix
        -> Int -> Int -- ^ size of result 2
        -> Ptr a -- ^ result ptr 3, a matrix
        -> Int -> Int -- ^ size of result 3
        -> Ptr a -- ^ matrix operand
        -> Int -- ^ CFFI call status code, non-zero means error
      ) -- ^ op function
    -> Mat a -- ^ operand
    -> (Mat a, Mat a, Mat a)
unsafeDecomposeOp (r1, c1) (r2, c2) (r3, c3) f m = unsafePerformIO $ do
    m1 <- Mutable.zeros r1 c1
    m2 <- Mutable.zeros r2 c2
    m3 <- Mutable.zeros r3 c3
    Mutable.unsafeWith m1 $ \vect1 _ _ ->
        Mutable.unsafeWith m2 $ \vect2 _ _ ->
            Mutable.unsafeWith m3 $ \vect3 _ _ ->
                unsafeWith m $ \xs _ _ -> I.call $ f vect1 r1 c1 vect2 r2 c2 vect3 r3 c3 xs
    m1' <- unsafeFreeze m1
    m2' <- unsafeFreeze m2
    m3' <- unsafeFreeze m3
    return (m1', m2', m3')

{-# INLINE unsafeDecomposeOp #-}

type instance Element (Mat a) = a

instance PrimType a => IsList (Mat a) where
    type Item (Mat a) = a
    fromList xs = M 1 (fromCount $ length xs) (fromList xs)
    fromListN n xs = M 1 n (fromListN n xs)
    toList = toList . vect

instance PrimType a => Collection (Mat a) where
    null M{..} = row == 0 && column == 0
    {-# INLINABLE null #-}
    length M{..} = CountOf (row * column)
    {-# INLINABLE length #-}
    elem x = elem x . vect
    {-# INLINABLE elem #-}
    minimum = minimum . nonEmpty_ . vect . getNonEmpty
    {-# INLINABLE minimum #-}
    maximum = maximum . nonEmpty_ . vect . getNonEmpty
    {-# INLINABLE maximum #-}
    all predicate = all predicate . vect
    {-# INLINABLE all #-}
    any predicate = any predicate . vect
    {-# INLINABLE any #-}

instance PrimType a => MutableCollection (Mutable.MMat a) where
    type MutableFreezed (Mutable.MMat a) = Mat a
    type MutableKey (Mutable.MMat a) = (Int, Int)
    type MutableValue (Mutable.MMat a) = a
    unsafeThaw (M r c v) = Mutable.M r c <$> unsafeThaw v
    unsafeFreeze (Mutable.M r c v) = M r c <$> unsafeFreeze v
    thaw (M r c v) = Mutable.M r c <$> thaw v
    freeze (Mutable.M r c v) = M r c <$> freeze v
    mutNew (CountOf c) = Mutable.M c 1 <$> mutNew (CountOf c)
    mutUnsafeWrite (Mutable.M _ c v) (i, j) = mutUnsafeWrite v (Offset (i * c + j))
    mutWrite (Mutable.M _ c v) (i, j) = mutWrite v (Offset (i * c + j))
    mutUnsafeRead (Mutable.M _ c v) (i, j) = mutUnsafeRead v (Offset (i * c + j))
    mutRead (Mutable.M _ c v) (i, j) = mutRead v (Offset (i * c + j))
