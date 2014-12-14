template <typename T>
void DMImage::GetData(std::vector<T> &destination, int t, int l, int b, int r, int front, int back)
{
	if (destination.size() != (r - l) * (b - t))
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
void DMImage::SetRealData(std::vector<T> &source)
{
	if (source.size() > width * height)
		throw std::invalid_argument("Array sizes must match");

	if (Image.IsDataTypeFloat())
		SetRealLockerData<float, T>(source);
	else if (Image.IsDataTypeInteger())
		SetRealLockerData<int, T>(source);
	else if (Image.IsDataTypeSignedInteger())
		SetRealLockerData<int, T>(source);
	else if (Image.IsDataTypeUnsignedInteger())
		SetRealLockerData<uint, T>(source);
	else
		throw std::invalid_argument("Data type not supported");
}

template <typename T>
void DMImage::SetComplexData(std::vector<T> &source, ShowComplex doComplex)
{

	if (source.size() > width * height){
		throw std::invalid_argument("Array sizes must match");
	}

	if (Image.IsDataTypeFloat())
		SetComplexLockerData<float, T>(source, doComplex);
	else if (Image.IsDataTypeInteger())
		SetComplexLockerData<int, T>(source, doComplex);
	else if (Image.IsDataTypeSignedInteger())
		SetComplexLockerData<int, T>(source, doComplex);
	else if (Image.IsDataTypeUnsignedInteger())
		SetComplexLockerData<uint, T>(source, doComplex);
	else
		throw std::invalid_argument("Data type not supported");
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
void DMImage::SetRealLockerData(std::vector<T> &source)
{

	Gatan::PlugIn::ImageDataLocker iLocker = Gatan::PlugIn::ImageDataLocker(Image);
	U* idata = static_cast<U*>(iLocker.get());
	for (int j = 0; j < (source.size()); j++)
		idata[j] = static_cast<U>(source[j]);
	iLocker.~ImageDataLocker();
	Image.DataChanged();
}

template <typename U, typename T>
void DMImage::SetComplexLockerData(std::vector<T> &source, ShowComplex doComplex)
{
	Gatan::PlugIn::ImageDataLocker iLocker = Gatan::PlugIn::ImageDataLocker(Image);
	U* idata = static_cast<U*>(iLocker.get());
	for (int j = 0; j < (source.size()); j++)
	{
		if (doComplex == SHOW_REAL)
			idata[j] = static_cast<U>(source[j].real());
		else if (doComplex == SHOW_IMAG)
			idata[j] = static_cast<U>(source[j].imag());
		else if (doComplex == SHOW_ABS)
			idata[j] = static_cast<U>(std::abs(source[j]));
	}
	iLocker.~ImageDataLocker();
	Image.DataChanged();
}