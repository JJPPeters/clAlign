template <typename T>
void DMImage::GetData(std::vector<T> &destination, int t, int l, int b, int r, int front, int back)
{
	if (destination.size() != (r - l) * (b - t) * (back - front))
		throw std::invalid_argument("Array sizes must match");

	if (Image.IsDataTypeFloat())
		GetLockerData<float, T>(destination, t, l, b, r, front, back);
	else if (Image.IsDataTypeInteger())
		GetLockerData<int, T>(destination, t, l, b, r, front, back);
	else if (Image.IsDataTypeSignedInteger())
		GetLockerData<int, T>(destination, t, l, b, r, front, back);
	else if (Image.IsDataTypeUnsignedInteger())
		GetLockerData<uint, T>(destination, t, l, b, r, front, back);
	else
		throw std::invalid_argument("Data type not supported");
}

template <typename T>
void DMImage::SetRealData(const std::vector<T> &source, int offset)
{
	if (source.size() > (width * height * depth) - offset)
		throw std::invalid_argument("Array sizes must match");

	if (Image.IsDataTypeFloat())
		SetRealLockerData<float, T>(source, offset);
	else if (Image.IsDataTypeInteger())
		SetRealLockerData<int, T>(source, offset);
	else if (Image.IsDataTypeSignedInteger())
		SetRealLockerData<int, T>(source, offset);
	else if (Image.IsDataTypeUnsignedInteger())
		SetRealLockerData<uint, T>(source, offset);
	else
		throw std::invalid_argument("Data type not supported");
}

template <typename T>
void DMImage::SetRealData(const std::vector<T> &source)
{
	SetRealData(source, 0);
}

template <typename T>
void DMImage::SetComplexData(const std::vector<T> &source, int offset, ShowComplex doComplex)
{
	if (source.size() > (width * height * depth) - offset){
		throw std::invalid_argument("Array sizes must match");
	}

	if (Image.IsDataTypeFloat())
		SetComplexLockerData<float, T>(source, offset, doComplex);
	else if (Image.IsDataTypeInteger())
		SetComplexLockerData<int, T>(source, offset, doComplex);
	else if (Image.IsDataTypeSignedInteger())
		SetComplexLockerData<int, T>(source, offset, doComplex);
	else if (Image.IsDataTypeUnsignedInteger())
		SetComplexLockerData<uint, T>(source, offset, doComplex);
	else
		throw std::invalid_argument("Data type not supported");
}

template <typename T>
void DMImage::SetComplexData(const std::vector<T> &source, ShowComplex doComplex)
{
	SetComplexData(source, 0, doComplex);
}

template <typename U, typename T>
void DMImage::GetLockerData(std::vector<T> &destination, int t, int l, int b, int r, int front, int back)
{
	Gatan::PlugIn::ImageDataLocker iLocker = Gatan::PlugIn::ImageDataLocker(Image);
	U* idata = static_cast<U*>(iLocker.get());
	for (int j = t; j < b; j++)
	for (int i = l; i < r; i++)
	for (int k = front; k < back; k++)
		destination[(i - t) + (j - t)*(r - l) + (k - front)*(r - l)*(b - t)] = static_cast<T>(idata[i + j * width + k * width * height]);
	iLocker.~ImageDataLocker();
}

template <typename U, typename T>
void DMImage::SetRealLockerData(const std::vector<T> &source, int offset)
{

	Gatan::PlugIn::ImageDataLocker iLocker = Gatan::PlugIn::ImageDataLocker(Image);
	U* idata = static_cast<U*>(iLocker.get());
	for (int j = 0; j < (source.size()); j++)
		idata[j+offset] = static_cast<U>(source[j]);
	iLocker.~ImageDataLocker();
	Image.DataChanged();
}

template <typename U, typename T>
void DMImage::SetComplexLockerData(const std::vector<T> &source, int offset, ShowComplex doComplex)
{
	Gatan::PlugIn::ImageDataLocker iLocker = Gatan::PlugIn::ImageDataLocker(Image);
	U* idata = static_cast<U*>(iLocker.get());
	for (int j = 0; j < (source.size()); j++)
	{
		if (doComplex == SHOW_REAL)
			idata[j+offset] = static_cast<U>(source[j].real());
		else if (doComplex == SHOW_IMAG)
			idata[j+offset] = static_cast<U>(source[j].imag());
		else if (doComplex == SHOW_ABS)
			idata[j+offset] = static_cast<U>(std::abs(source[j]));
	}
	iLocker.~ImageDataLocker();
	Image.DataChanged();
}