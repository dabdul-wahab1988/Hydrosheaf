# 🚀 Hydrosheaf Web Application - Quick Start Guide

**Critical Issue: FIXED!** ✅ The backend now uses real Hydrosheaf analysis instead of mock data.

## 30-Second Setup

```bash
# 1. Install Hydrosheaf in backend
cd web/backend
pip install -r requirements.txt
pip install ../../  # ← CRITICAL: Install Hydrosheaf core

# 2. Start backend (Terminal 1)
uvicorn app.main:app --reload

# 3. Start frontend (Terminal 2)
cd ../frontend
npm install
npm run dev

# 4. Test (Terminal 3)
cd ../backend
python test_integration.py
```

## Open in Browser

- **Frontend**: http://localhost:5173
- **Backend API**: http://localhost:8000/api/docs

## Try It Out

1. Go to http://localhost:5173
2. Click **"Analysis"** in sidebar
3. Select **"Demo Groundwater Samples"** from dropdown
4. Click **"Start Analysis"**
5. Wait 5-10 seconds
6. See **REAL Hydrosheaf results**! 🎉

## Verify Real Analysis

Check the results metadata:
```json
{
  "metadata": {
    "analysis_engine": "Hydrosheaf Core",
    "mock_data": false  ← Should be FALSE!
  }
}
```

## What Works Now

✅ Real transport modeling (evaporation vs mixing)
✅ Real LASSO reaction fitting (sparse selection)
✅ Real network flow inference (Bayesian)
✅ Automatic data transformation (mg/L → mmol/L)
✅ Feature auto-disable (missing data handling)

## Troubleshooting

### "Hydrosheaf not available"

```bash
cd web/backend
pip install ../../
python -c "import hydrosheaf"  # Should not error
```

### Backend won't start

```bash
# Check Python version
python --version  # Need >= 3.8

# Reinstall dependencies
pip install -r requirements.txt
```

### Frontend won't start

```bash
# Check Node version
node --version  # Need >= 18.x

# Reinstall dependencies
rm -rf node_modules package-lock.json
npm install
```

## Full Documentation

- **Setup Details**: `web/backend/README.md`
- **What Was Fixed**: `INTEGRATION_FIXED.md`
- **Full Analysis**: `WEB_INTEGRATION_ANALYSIS.md`
- **Integration Guide**: `web/backend/INTEGRATION_GUIDE.md`

## Need Help?

Run the test suite:
```bash
cd web/backend
python test_integration.py
```

If all tests pass → Everything is working! 🎉

If tests fail → Check error messages and see troubleshooting above.

---

**Ready to do real hydrogeochemistry!** 💧🔬
